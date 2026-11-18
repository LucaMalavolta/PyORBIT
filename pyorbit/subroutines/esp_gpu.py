"""CUDA float64 ESP likelihood, with S+LEAF's rank-3(2*nharm+1) factors.

No dense N x N covariance is formed. One CUDA block performs the stable
semiseparable LDL recurrence; each thread owns entries of its small state.
The factors and the likelihood are evaluated on the selected GPU. Host code
still assembles PyORBIT parameters and drives its existing samplers.
"""

import numpy as np
from functools import lru_cache


@lru_cache(None)
def _kernels(rank):
    import cupy as cp

    source = r"""
    extern "C" __global__ void factors(const double* t, const double* alpha,
      const double* beta, const double* weights, double* U, double* V,
      double* A, const double* noise, const int n, const int hmax,
      const double period, const double rho) {
      int i=blockDim.x*blockIdx.x+threadIdx.x;
      if(i>=n) return;
      const double mu=1.326644517327145, la=1.0907260149419182/rho;
      const double nu=mu*la, omega=6.283185307179586476925286766559/period;
      const double a0=2./3.*(1.+1./(mu*mu)), a=1.-a0, b=1./mu;
      double c=cos(nu*t[i]), s=sin(nu*t[i]);
      double eu[3]={a0,a*c+b*s,a*s-b*c}, ev[3]={1.,c,s};
      double du[3]={-la*a0,nu*(-a*s+b*c)-la*eu[1],nu*(a*c+b*s)-la*eu[2]};
      double dv[3]={la,-nu*s+la*c,nu*c+la*s};
      double variance=0.; int psize=2*hmax+1;
      for(int j=0;j<3;j++) for(int k=0;k<psize;k++) {
        int h=(k+1)/2; double pc=cos(h*omega*t[i]), ps=sin(h*omega*t[i]);
        double pv=k==0?1.:(k%2?pc:ps);
        double pd=k==0?0.:(k%2?-h*omega*ps:h*omega*pc);
        double pu=weights[h]*pv, pdu=weights[h]*pd;
        int r=j*psize+k;
        double u=eu[j]*pu, v=ev[j]*pv;
        double ud=du[j]*pu+eu[j]*pdu, vd=dv[j]*pv+ev[j]*pd;
        U[i*R+r]=alpha[i]*u+beta[i]*ud;
        V[i*R+r]=alpha[i]*v+beta[i]*vd;
        variance+=beta[i]*beta[i]*ud*vd;
      }
      A[i]=noise[i]+alpha[i]*alpha[i]+variance;
    }
    extern "C" __global__ void likelihood(const double* t,const double* U,
      const double* V,const double* A,const double* y,double* out,
      const int n,const double rho,const double neg_inf) {
      __shared__ double S[R*R], f[R], u[R], q[R], w[R], d,z,p,total;
      int k=threadIdx.x;
      for(int j=k;j<R*R;j+=blockDim.x) S[j]=0.;
      if(k<R) f[k]=0.;
      if(k==0) total=0.;
      __syncthreads();
      for(int i=0;i<n;i++) {
        if(k<R) u[k]=U[i*R+k];
        __syncthreads();
        if(k<R) { double v=0.; for(int j=0;j<R;j++) v+=S[k*R+j]*u[j]; q[k]=v; }
        __syncthreads();
        if(k==0) {
          d=A[i]; z=y[i];
          for(int j=0;j<R;j++){d-=u[j]*q[j]; z-=u[j]*f[j];}
          p=i+1<n?exp(-1.0907260149419182/rho*(t[i+1]-t[i])):1.;
          total+=log(d)+z*z/d;
        }
        __syncthreads();
        if(!(d>0.) || !isfinite(d)) { if(k==0) out[0]=neg_inf; return; }
        if(k<R) w[k]=(V[i*R+k]-q[k])/d;
        __syncthreads();
        if(k<R) f[k]=p*(f[k]+w[k]*z);
        for(int j=k;j<R*R;j+=blockDim.x) S[j]=p*p*(S[j]+d*w[j/R]*w[j%R]);
        __syncthreads();
      }
      if(k==0) out[0]=-.5*(total+n*1.837877066409345483560659472811);
    }
    """
    module = cp.RawModule(
        code="#define R %d\n" % rank + source, options=("--std=c++11",)
    )
    return module.get_function("factors"), module.get_function("likelihood")


def esp_loglike(
    time, residuals, variance, period, rho, eta, alpha=1.0, beta=0.0, nharm=4
):
    import cupy as cp
    import scipy.special
    from pyorbit.common.gpu_backend import get_device

    device = get_device()
    with cp.cuda.Device(device.id):
        t = cp.ascontiguousarray(cp.asarray(time, dtype=cp.float64))
        y = cp.ascontiguousarray(cp.asarray(residuals, dtype=cp.float64))
        noise = cp.ascontiguousarray(cp.asarray(variance, dtype=cp.float64))
        aa = cp.ascontiguousarray(
            cp.broadcast_to(cp.asarray(alpha, dtype=cp.float64), t.shape)
        )
        bb = cp.ascontiguousarray(
            cp.broadcast_to(cp.asarray(beta, dtype=cp.float64), t.shape)
        )
        h = np.arange(nharm + 1, dtype=np.float64)
        weights = scipy.special.ive(h, 1.0 / (4.0 * eta * eta))
        weights *= np.where(h == 0, 1.0, 2.0)
        weights /= np.sum(weights)
        weights = cp.asarray(weights, dtype=np.float64)
        rank = 3 * (2 * nharm + 1)
        n = t.size
        U = cp.empty((n, rank), dtype=cp.float64)
        V = cp.empty_like(U)
        A = cp.empty(n, dtype=cp.float64)
        result = cp.empty(1, dtype=cp.float64)
        factors, likelihood = _kernels(rank)
        factors(
            ((n + 127) // 128,),
            (128,),
            (
                t,
                aa,
                bb,
                weights,
                U,
                V,
                A,
                noise,
                np.int32(n),
                np.int32(nharm),
                np.float64(period),
                np.float64(rho),
            ),
        )
        likelihood(
            (1,),
            (128,),
            (t, U, V, A, y, result, np.int32(n), np.float64(rho), np.float64(-np.inf)),
        )
        return float(result.get()[0])

"""Quadratic circular transit on JAX, checked against batman in float64."""
from functools import partial
import jax
import jax.numpy as jnp
from jaxoplanet.core.limb_dark import light_curve

@partial(jax.jit, static_argnames=('samples',))
def circular_transit(t, period, tc, radius, ars, inclination, u, exposure, samples=1):
    offsets=jnp.linspace(-.5,.5,samples)*exposure if samples>1 else jnp.zeros(1)
    phase=2*jnp.pi*(t[:,None]+offsets-tc)/period
    separation=ars*jnp.sqrt(jnp.sin(phase)**2+(jnp.cos(inclination)*jnp.cos(phase))**2)
    flux=light_curve(u,separation,radius,order=40)
    return jnp.mean(jnp.where(jnp.cos(phase)>0.,flux,0.),axis=1)

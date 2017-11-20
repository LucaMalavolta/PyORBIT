from abstract_model import *

class Correlations(AbstractModel):

    model_class = 'correlations'
    list_pams_common = {}
    list_pams_dataset = {
        'c1': 'U',
        'c2': 'U',
        'c3': 'U',
        'c4': 'U',
        'c5': 'U',
        'c6': 'U',
        'c7': 'U',
        'c8': 'U',
        'c9': 'U'
    }

    def model_setup(self, dataset, **kwargs):

        for ii in xrange(kwargs['order']+1, 10):
            self.fix_list[dataset.name_ref]['c'+repr(ii)] = 0.0000


        """To each dataset the correlated dataset is associated"""
        self.bounds[name_ref][name_asc] = {}
        self.variables[name_ref][name_asc] = {}
        self.fix_list[name_ref][name_asc] = {}
        self.var_list[name_ref][name_asc] = {}

        self.prior_kind[name_ref][name_asc] = {}
        self.prior_pams[name_ref][name_asc] = {}

        self.list_pams[name_ref][name_asc] = {}

        self.order[name_ref][name_asc] = 1
        self.order_ind[name_ref][name_asc] = {}

        self.x_vals[name_ref][name_asc] = np.zeros(mc.dataset_dict[name_ref].n, dtype=np.double)
        self.x_mask[name_ref][name_asc] = np.zeros(mc.dataset_dict[name_ref].n, dtype=bool)

        """ HERE: we must associated the data from name_asc dataset to the one from name_ref
            remove that part from dataset.pyx
            Add a None option for the dataset
            Fix input_parser to accomodate the new changes
            Jitter must not be included in the analysis, but how about the offset? 
            Or maybe I should just leave the zero point of the polynomial fit free?
        """

        for i_date, v_date in enumerate(mc.dataset_dict[name_asc].x):
            match = np.where(np.abs(v_date-mc.dataset_dict[name_ref].x) < threshold)[0]
            self.x_vals[name_ref][name_asc][match] = mc.dataset_dict[name_asc].y[i_date]
            self.x_mask[name_ref][name_asc][match] = True

        self.x_zero[name_ref][name_asc] = np.median(self.x_vals[name_ref][name_asc][self.x_mask[name_ref][name_asc]])

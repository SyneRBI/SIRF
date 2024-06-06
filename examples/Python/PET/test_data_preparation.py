import os
from sirf.Utilities import examples_data_path

import importlib
pet = importlib.import_module('sirf.STIR')

data_path = examples_data_path('PET') + '/mMR'
print('Finding files in %s...' % data_path)

#save_path = "~/tmp/data"
#print('Saving files in %s...' % save_path)

f_listmode = os.path.join(data_path, "list.l.hdr");
f_template = os.path.join(data_path, "mMR_template_span11_small.hs");
f_attn = os.path.join(data_path, "mu_map.hv");
f_norm = os.path.join(data_path, "norm.n.hdr");

# engine's messages go to files, except error messages, which go to stdout
_ = pet.MessageRedirector('info.txt', 'warn.txt')

# select acquisition data storage scheme
pet.AcquisitionData.set_storage_scheme('memory')

# read acquisition data template
acq_data_template = pet.AcquisitionData(f_template)

listmode_data = pet.ListmodeData(f_listmode)

# create listmode-to-sinograms converter object
lm2sino = pet.ListmodeToSinograms()
lm_data = pet.ListmodeData(f_listmode)
prompts, randoms = lm2sino.prompts_and_randoms_from_listmode(lm_data, 0, 10, acq_data_template)
print('data shape: %s' % repr(prompts.shape))
print('prompts norm: %f' % prompts.norm())
print('randoms norm: %f' % randoms.norm())
#prompts.write(os.path.join(save_path, 'prompts.hs'))
#randoms.write(os.path.join(save_path, 'randoms.hs'))

attn_image = pet.ImageData(f_attn)
attn, acf, iacf = pet.AcquisitionSensitivityModel.compute_attenuation_factors(prompts, attn_image)
print('norm of the attenuation correction factor: %f' % acf.norm())
print('norm of the inverse of the attenuation correction factor: %f' % iacf.norm())
#acf.write(os.path.join(save_path, 'acf.hs'))
#iacf.write(os.path.join(save_path, 'iacf.hs'))
#exit()

asm = pet.AcquisitionSensitivityModel(f_norm)
se = pet.ScatterEstimator()
se.set_input(prompts)
se.set_attenuation_image(attn_image)
se.set_randoms(randoms)
se.set_asm(asm)
se.set_attenuation_correction_factors(iacf)
se.set_num_iterations(4)
se.set_OSEM_num_subsets(7)
se.set_output_prefix("scatter")
se.set_up()
se.process()
scatter = se.get_output()
print('norm of the scatter estimate: %f' % scatter.norm())

multfact = acf.clone()
asm.set_up(acf)
asm.unnormalise(multfact)
print(multfact.norm())

background = randoms + scatter
print('norm of the backgrount term: %f' % background.norm())

asm_mf = pet.AcquisitionSensitivityModel(multfact)
asm_mf.set_up(background)
asm_mf.normalise(background)
print('norm of the additive term: %f' % background.norm())


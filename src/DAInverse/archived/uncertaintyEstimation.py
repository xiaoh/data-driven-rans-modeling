
# Initiate an instance of ReynoldsStressRF

ReynoldsStressRF reynolds()
# give dir of OF,
# give number of realizations (Nens = 5);


# Initiasize an instance of KLReducedModel
KLReducedModel deltaXi()
KLReducedModel deltaEta()
## give Nens

# Initiasize an instance of enKF
enKF enkfFoam();
# Give Nens
# case folders

enkfFoam.initialize();
# propagate tau ensembles to velocity

# select another set of omega_\xi, omega_eta
# (in practice, we get this from filtered states)

# update omega_\xi, omega_eta here



deltaXi.recomputeField();
deltaEta.recomputeField();

# update deltaXi and deltaEta realizations
reynolds.modifyError(deltaXi.field(), deltaEta.field())

reynolds.perturbTau()
# write perturbed fields to modified cases

enkfFoam.propgate();
# To obtain new enssemble of velocities.



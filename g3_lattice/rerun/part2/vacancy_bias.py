"""Central vacancy-pair feature geometry; no coefficients are guessed here.
F_i=sum_s beta_s sum_{vacant j in shell s} Rhat_ij, beta in eV/A.
All physical periodic images are included, unlike a single minimum-image tie.
The source/fit occupancy is a bulk reference pattern; an open displacement
boundary is NOT reinterpreted as additional vacancies.
"""
import itertools
from math import sqrt
import numpy as np


SHELL_CART2_SQUARED={1:2,2:4,3:6}


def central_vacancy_features(occupied,CELL,shells=(1,2)):
    occupied=np.asarray(occupied)
    if occupied.ndim!=3 or occupied.dtype!=np.bool_:
        raise ValueError('A three-dimensional reference occupation mask is required')
    basis2=2*np.asarray(CELL)
    if not np.array_equal(basis2,np.rint(basis2)):
        raise ValueError('Expected the delivered half-integer FCC basis')
    targets=[SHELL_CART2_SQUARED[s] for s in shells]
    features=np.zeros((len(shells),3)+occupied.shape)
    counts=[0]*len(shells)
    for offset in itertools.product(range(-3,4),repeat=3):
        cart2=basis2@offset;squared=int(np.dot(cart2,cart2))
        if squared not in targets:continue
        shell=targets.index(squared)
        direction=cart2/sqrt(squared)
        vacant_neighbor=np.roll(~occupied,tuple(-v for v in offset),axis=(0,1,2))
        features[shell]+=direction[:,None,None,None]*vacant_neighbor
        counts[shell]+=1
    expected={1:12,2:6,3:24}
    if counts!=[expected[s] for s in shells]:
        raise ArithmeticError('Physical FCC shell enumeration is incomplete')
    return features,{'shells':list(shells),'physical_image_neighbors_per_shell':counts,
        'units':'unit-vector sums; beta_s in eV/A per physical vacant-neighbor image',
        'boundary':'periodic reference occupation, independent of open/clamped displacement exterior',
        'half_cell_ties':'Both opposite physical images retained; odd force contributions may cancel exactly',
        'rank_scope':'Zero/unidentified columns must not be reported as measured zero coefficients'}


def predict_reference_bias(occupied,CELL,model,second_shell_multiplier=0.):
    """Constant reference-force term; no spring derivative is added again."""
    if model['schema']!='central-vacancy-bias-v1':
        raise ValueError('Unknown fitted bias model')
    if second_shell_multiplier not in (0.,.5,-.5):
        raise ValueError('Only the prescribed0,+C1/2,-C1/2 scenarios are in scope')
    features,contract=central_vacancy_features(occupied,CELL,(1,2,3))
    coefficients=model['shell_coefficients_eV_A']
    beta=np.array([coefficients[str(s)] for s in (1,2,3)],float)
    if not np.isfinite(beta).all() or beta[1]!=0:
        raise ValueError('Base model must keep director819 assumed C2=0 explicit')
    shell2_norm=float(np.linalg.norm(features[1,:,occupied]))
    if second_shell_multiplier and shell2_norm==0:
        raise ValueError('Sensitivity cannot use a2a-periodic occupancy with zero shell2 feature')
    beta[1]=second_shell_multiplier*beta[0]
    force=np.einsum('s,sb...->b...',beta,features)*occupied
    metadata=dict(kind='fitted_central_vacancy_law',authority='812/819',
        coefficients_used_eV_A={str(s):float(v) for s,v in zip((1,2,3),beta)},
        second_shell_multiplier=second_shell_multiplier,second_shell_feature_L2=shell2_norm,
        second_shell_scope='Base0 is an assumption; +/-halfC1 are scenarios, not fitted bounds',
        force_origin='Declared reference sites, not relaxed/F0-over-K shifted',
        force_definition='Constant F0; directions are reference lines, not recomputed from instantaneous u',
        boundary=contract['boundary'],accuracy_scope=model['accuracy_scope'])
    return force,metadata

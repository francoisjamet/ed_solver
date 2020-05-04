#!/bin/env python2.7
import numpy as np
import h5py
import subprocess as sp
def create_file(beta,nomg,U,U1,e1,e2,t) :
    W = np.array([ 1j*(2*i+1)*np.pi/beta for i in range(nomg) ])
    D1 = 0.2/W
    D2 = 0.3/W
    D = np.zeros((nomg,2,2), dtype= complex)
    D[:,0,0] = D1
    D[:,1,1] = D2
    f = open('delta1.inp','w')
    ff = open('delta2.inp','w')
    for iomg in range(nomg) :
        s = str(W[iomg].imag) +' '
        for i in range(2):
            for j in range(2):
                d = D[iomg,i,j]
                s += ' {} {} '.format(d.real,d.imag)
        s+='\n'
        f.write(s)
        ff.write(s)
    f.close()
    ff.close()

    f = open('PARAMS','w')
    s = """2
{e1} {t}
{t} {e2}
{e1} {t}
{t} {e2}
{U}  {U1}
{U1} {U}
1000
1000
1000
F
1
F
0.0000000000000000
0.0000000000000000
0.
{beta}
0.0    !8 UU
1.0000000000000001E-002
3.4000000000000000E-001
3.2999999999999999E-001
   7.0999999999999997E-001
          10
  0.000000000000000  !0.7 ! JJ
 F
 F
 F
 F
0.0
    """.format(U=U,U1=U1,e1=e1,e2=e2,t=t,beta=beta)
    f.write(s)
    f.close()



    f = open('cutoff','w')
    f.write("""1
2
1d-4
1d-6""")
    f.close()


    sp.call('mkdir -p ED',shell=True)
    f = open('ED/ED.in','w')
    f.write("""


 PAIRING_IMP_TO_BATH=.false.
 track_sectors=.true.
 fast_fit=.true.

 first_iter_use_edinput=.false.

 start_para=.false.
 FLAG_GUP_IS_GDN=.false.

 force_nupdn_basis=.true.
 force_sz_basis=.false.
 force_no_pairing=.true.

 lambda_sym_fit=0.p
 fit_shift=0.0d0
 weight_expo=2

 FLAG_MPI_GREENS=0
 window_hybrid=0
 window_hybrid2=0
 window_weight=0.
 search_step=0.000001

 cutoff_hamilt_param=0.000000010
 dist_max=1.d-12
 tolerance=1.d-13
 FIT_METH=CIVELLI

 flag_introduce_only_noise_in_minimization=.false.
 flag_introduce_noise_in_minimization=.false.


flag_idelta_two_scales_ed=0
nsec=-1
nsec0=0
which_lanczos=NORMAL
Neigen=4
Block_size=0
ncpt_approx=0
cpt_upper_bound=1.00000000
cpt_lagrange=0.00000001
iwindow=1
fmos_iter=10
fmos_mix=0.40000000
fmos=.false.
fmos_fluc=.false.
fmos_hub1=.false.
ncpt_flag_two_step_fit=.false.
Niter_search_max_0=800000
FLAG_ALL_GREEN_FUNC_COMPUTED=.true.
fit_weight_power=1.00000000
force_para_state=.true.
FLAG_GUP_IS_GDN=.true.
fit_nw=200
min_all_bath_param=3

FLAG_DUMP_INFO_FOR_GAMMA_VERTEX=.false.
Nitermax=1000
Nitergreenmax=400

Neigen=5

which_lanczos=FULL_ED
FLAG_FULL_ED_GREEN=.false.
dEmax0=10
Block_size=0
tolerance=1.d-13
    """)
    f.close()
    sp.call('cd ED && generate_edcorrel  2 2 && mv ed.correl  ed_correl1 ',shell=True)

# def run_triqs () :
#     # This part run triqs with the same model for testing.
#     import pytriqs.utility.mpi as mpi
#     from pytriqs.gf import *
#     from pytriqs.operators import *
#     from pytriqs.applications.impurity_solvers.cthyb import Solver
#     from pytriqs.gf.tools import *
#     from pytriqs.archive import HDFArchive

#     Eimp = np.array([[e1, t],
#                      [t,e2]])
#     I = np.identity(2)
#     G0 = 0 * D
#     for iw in range(nomg) :
#         G0[iw] = np.linalg.inv((W[iw] * I -D[iw] - Eimp))


#         # Construct the impurity solver with the inverse temperature
#         # and the structure of the Green's functions
#     S = Solver(beta = beta, gf_struct = {'up':[0,1], 'down':[0,1]},n_iw=nomg)

#     # Initialize the non-interacting Green's function S.G0_iw

#     for spin, g0 in S.G0_iw:
#         g0.data[nomg:,:,:] = G0[:,:,:]
#         g0.data[:nomg,:,:] = np.conjugate(G0[::-1,:,:])

#     params = {}
#     params['n_cycles'] = 30000                # Number of QMC cycles
#     params['length_cycle'] = 2000                # Length of one cycle
#     params['n_warmup_cycles'] = 10000           # Warmup cycles
#     params['measure_G2_iw'] = True
#     params["measure_G2_n_fermionic"]=10


#     # Run the solver. The result will be stored in S.G_tau.
#     H  = U * n('up',0) * n('down',0) +  U * n('up',1) * n('down',1)
#     H += U1 * n('up',0) * n('down',1) +  U1 * n('up',0) * n('down',1)
#     H += U1 * n('up',0) * n('up',1) +  U1 * n('up',0) * n('up',1)
#     S.solve(h_int = H, **params)

#     if mpi.is_master_node():
#         with HDFArchive('imp.h5','w') as A:
#             A['S']=S
#             A['G']=S.G_iw

#     for spin, g in S.G_iw:
#         write_gf_to_txt(g)


def run(beta) :
    sp.call('dmft_solver > out1',shell=True)
    sp.call('cp sig1.inp sig',shell=True)
    sp.call('cp g1.inp g',shell=True)
    sp.call('dmft_solver FLAG_DUMP_INFO_FOR_GAMMA_VERTEX=.true. dEmax0=10000 Neigen=500 Nitergreenmax=1 >out2',shell=True)
    sp.call('omega_path 10',shell=True)
    sp.call('env OMP_NUM_THREADS=3  dmft_chiloc_2 >chiloc.out',shell=True)

    sp.call('read_g2 2 {}'.format(beta),shell=True)
    # cs = np.array(h5py.File('chiS.h5','r')['chis'])
    # cc = np.array(h5py.File('chiS.h5','r')['chic'])
    # g = np.loadtxt('g')
    # csb = np.array(h5py.File('../bm/chiS_bm.h5','r')['chis'])
    # ccb = np.array(h5py.File('../bm/chiS_bm.h5','r')['chic'])
    # gb = np.loadtxt('../bm/g_bm')
    # if np.allclose(cs,csb) and np.allclose(cs,csb)  :
    #     print('test passed')
    # else :
    #     print('test failed')

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description='create and run test for ed. Make sure you have the directory bin in your file')
    p.add_argument('--run',action='store_true')
    args = p.parse_args()
    run_test = args.run



    U = 10
    beta = 50
    e1=  -U / 2.
    e2 = -U / 2.
    t = 0
    U1 = 2
    nomg = 1000

    create_file(beta,nomg,U,U1,e1,e2,t)

    if run_test :
        run(beta)

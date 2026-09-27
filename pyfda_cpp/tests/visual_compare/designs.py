# shared battery: (name, pyfda fil[0] updates, pyfda_cpp spec)
import math
def a_pb_lin(db, fir=False):  # pyfda linear pass band spec
    return (10**(db/20)-1)/(10**(db/20)+1) if fir else 1 - 10**(-db/20)
def a_sb_lin(db): return 10**(-db/20)
D = [
 ('d1_ellip_lp', dict(ft='IIR', fc='Ellip', rt='LP', fo='man', N=4, f_c=0.1, f_pb=0.1, f_sb=0.2, A_PB=2, A_SB=60)),
 ('d2_cheby1_hp_min', dict(ft='IIR', fc='Cheby1', rt='HP', fo='min', N=4, f_pb=0.3, f_sb=0.2, A_PB=1, A_SB=40)),
 ('d3_butter_bp', dict(ft='IIR', fc='Butter', rt='BP', fo='man', N=4, f_c=0.2, f_c2=0.3, f_pb=0.2, f_pb2=0.3, f_sb=0.15, f_sb2=0.35, A_PB=3, A_SB=40)),
 ('d4_equiripple_lp', dict(ft='FIR', fc='Equiripple', rt='LP', fo='man', N=40, f_pb=0.1, f_sb=0.15, A_PB=1, A_SB=60)),
 ('d5_firwin_bs', dict(ft='FIR', fc='Firwin', rt='BS', fo='man', N=60, f_c=0.2, f_c2=0.3, f_pb=0.15, f_pb2=0.35, f_sb=0.2, f_sb2=0.3, A_PB=1, A_SB=60)),
 ('d6_cheby2_lp_min', dict(ft='IIR', fc='Cheby2', rt='LP', fo='min', N=4, f_pb=0.1, f_sb=0.15, A_PB=1, A_SB=60)),
]

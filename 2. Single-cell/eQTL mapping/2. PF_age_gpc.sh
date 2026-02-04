 #!/bin/bash

module load PEER

peertool -f Adv_fib_peer.csv -c Adv_fib_age_pc.csv -n 20 -i 500 -o ./Adv_fib_age_pc
peertool -f Alv_fib_peer.csv -c Alv_fib_age_pc.csv -n 20 -i 500 -o ./Alv_fib_age_pc
peertool -f Alv_mph_peer.csv -c Alv_mph_age_pc.csv -n 20 -i 500 -o ./Alv_mph_age_pc
peertool -f Alv_trans_peer.csv -c Alv_trans_age_pc.csv -n 20 -i 500 -o ./Alv_trans_age_pc
peertool -f AT1_peer.csv -c AT1_age_pc.csv -n 20 -i 500 -o ./AT1_age_pc
peertool -f AT2_peer.csv -c AT2_age_pc.csv -n 20 -i 500 -o ./AT2_age_pc

peertool -f Basal_peer.csv -c Basal_age_pc.csv -n 20 -i 500 -o ./Basal_age_pc
peertool -f Bcells_peer.csv -c Bcells_age_pc.csv -n 20 -i 500 -o ./Bcells_age_pc
peertool -f CD4_peer.csv -c CD4_age_pc.csv -n 20 -i 500 -o ./CD4_age_pc
peertool -f CD8_peer.csv -c CD8_age_pc.csv -n 20 -i 500 -o ./CD8_age_pc
peertool -f Cla_mono_peer.csv -c Cla_mono_age_pc.csv -n 20 -i 500 -o ./Cla_mono_age_pc
peertool -f Club_peer.csv -c Club_age_pc.csv -n 20 -i 500 -o ./Club_age_pc

peertool -f DC1_peer.csv -c DC1_age_pc.csv -n 20 -i 500 -o ./DC1_age_pc
peertool -f DC2_peer.csv -c DC2_age_pc.csv -n 20 -i 500 -o ./DC2_age_pc
peertool -f EC_aero_cap_peer.csv -c EC_aero_cap_age_pc.csv -n 20 -i 500 -o ./EC_aero_cap_age_pc
peertool -f EC_art_peer.csv -c EC_art_age_pc.csv -n 20 -i 500 -o ./EC_art_age_pc
peertool -f EC_gen_cap_peer.csv -c EC_gen_cap_age_pc.csv -n 20 -i 500 -o ./EC_gen_cap_age_pc
peertool -f EC_ven_pul_peer.csv -c EC_ven_pul_age_pc.csv -n 20 -i 500 -o ./EC_ven_pul_age_pc

peertool -f EC_ven_sys_peer.csv -c EC_ven_sys_age_pc.csv -n 20 -i 500 -o ./EC_ven_sys_age_pc
peertool -f Goblet_peer.csv -c Goblet_age_pc.csv -n 20 -i 500 -o ./Goblet_age_pc
peertool -f Int_mph_peri_peer.csv -c Int_mph_peri_age_pc.csv -n 20 -i 500 -o ./Int_mph_peri_age_pc
peertool -f Lym_EC_mat_peer.csv -c Lym_EC_mat_age_pc.csv -n 20 -i 500 -o ./Lym_EC_mat_age_pc
peertool -f Lym_EC_pro_peer.csv -c Lym_EC_pro_age_pc.csv -n 20 -i 500 -o ./Lym_EC_pro_age_pc
peertool -f Mast_peer.csv -c Mast_age_pc.csv -n 20 -i 500 -o ./Mast_age_pc

peertool -f Mig_DC_peer.csv -c Mig_DC_age_pc.csv -n 20 -i 500 -o ./Mig_DC_age_pc
peertool -f Mono_mph_peer.csv -c Mono_mph_age_pc.csv -n 20 -i 500 -o ./Mono_mph_age_pc
peertool -f Multiciliated_peer.csv -c Multiciliated_age_pc.csv -n 20 -i 500 -o ./Multiciliated_age_pc
peertool -f NK_peer.csv -c NK_age_pc.csv -n 20 -i 500 -o ./NK_age_pc
peertool -f Noncla_mono_peer.csv -c Noncla_mono_age_pc.csv -n 20 -i 500 -o ./Noncla_mono_age_pc
peertool -f Peri_fib_peer.csv -c Peri_fib_age_pc.csv -n 20 -i 500 -o ./Peri_fib_age_pc

peertool -f Sec_trans_peer.csv -c Sec_trans_age_pc.csv -n 20 -i 500 -o ./Sec_trans_age_pc
peertool -f SM_peer.csv -c SM_age_pc.csv -n 20 -i 500 -o ./SM_age_pc
peertool -f Sub_fib_peer.csv -c Sub_fib_age_pc.csv -n 20 -i 500 -o ./Sub_fib_age_pc



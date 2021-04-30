       integer,parameter :: pid_synchotron=1
       integer,parameter :: pid_nn_core=2
       integer,parameter :: pid_np_core=4
       integer,parameter :: pid_pp_core=8
       integer,parameter :: pid_nn_inner_crust=16
       integer,parameter :: pid_nn_core_super=32
       integer,parameter :: pid_np_core_super=64
       integer,parameter :: pid_pp_core_super=128
       integer,parameter :: pid_nn_inner_crust_super=256
       integer,parameter :: pid_PBF_s_p_core=512
       integer,parameter :: pid_PBF_s_n_core=1024
       integer,parameter :: pid_PBF_s_n_inner_crust=2048
       integer,parameter :: pid_PBF_pA_core=4096
       integer,parameter :: pid_PBF_pB_core=8192
       integer,parameter :: pid_PBF_pA_inner_crust=16384
       integer,parameter :: pid_PBF_pB_inner_crust=32768
       integer,parameter :: pid_mp_core=65536
       integer,parameter :: pid_ep_core=131072
       integer,parameter :: pid_eI_crust=262144
       integer,parameter :: pid_B1=524288
       integer,parameter :: pid_B2=1048576
       integer,parameter :: pid_B3=2097152

       
       real*8 :: gann
       real*8 :: gapp
       real*8 :: gaee
       real*8 :: gamm
       common/Couplings/gann,gapp,gaee,gamm

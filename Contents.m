% FDASRVF
%
% Files
%   align_fPCA                    - Group-wise function alignment and PCA Extractions
%   ampbox                        - ampbox Produce amplitude box plots
%   bootTB                        - Boostrapped tolerance bounds
%   calculateCentroid             - Calculate centroid of curve
%   create_basismatrix            - Create Basis Matrix
%   curve_to_q                    - Convert curve to Square-Root Velocity Function        
%   elastic_depth                 - compute elastic depth   
%   elastic_distance              - Calculates the two elastic distances between two functions
%   elastic_distance_curve        - Calculates the elastic distance between two curves
%   elastic_logistic              - elastic_logistic A class to provide SRVF logistic regression
%   elastic_lcpr_regression       - elastic logistic principal component regression
%   elastic_mlogistic             - elastic_mlogistic A class to provide SRVF multionomial logistic
%   elastic_mlpcr_regression      - elastic multinmoical principal component regression
%   elastic_regression            - elastic_regression A class to provide SRVF regression
%   elastic_pcr_regression        - elastic principal component regression
%   exp_map                       - Exponential Map
%   f_to_srvf                     - Convert function to Square-Root Velocity Function
%   fdacurve                      - fdacurve A class to provide registration of curves in R^n using SRVF
%   fdahpca                       - fdahpca A class to provide horizontal fPCA
%   fdahpns                       - fdahpns A class to pronivde horizontal fPNS
%   fdajpca                       - fdajpca A class to provide joint fPCA
%   fdakma                        - fdakma A class to provide a kmeans clustering and alignment
%   fdavpca                       - fdavpca A class to provide vertical fPCA
%   fdawarp                       - fdawarp A class to provide alignment of functional data using SRVF
%   gam_to_h                      - convert warping function to h space
%   gam_to_psi                    - convert warping fucntion to hilbert sphere
%   gam_to_v                      - convert warping function to tangent space on hilbert sphere
%   geodesic_sphere_Full          - Calculates geodesic on sphere
%   h_to_gam                      - convert h to warping function
%   inv_exp_map                   - Inverse Exponential Map
%   invertGamma                   - Invert Warping Function
%   L2norm                        - L2 Functional Norm
%   optimum_reparam               - Calculates Warping for two SRVFs
%   optimum_reparam_curve         - Calculates warping for two curve SRVFs
%   outlier_detection             - Outlier Detection
%   pairwise_align                - Align two functions
%   pairwise_align_bayes          - PAIRWISE_ALIGN Align two functions using Bayesian method
%   pairwise_align_curves         - registers to curves
%   Path_Plot                     - Plot geodesic path
%   pcaTB                         - functional principal component tolerance bound generation
%   phbox                         - phbox Construct phase box plots
%   psi_to_gam                    - Hilbert sphere to warping function
%   q_to_curve                    - Convert SRVF to curve
%   ReSampleCurve                 - Resample curve to have N points
%   rgam                          - Generate random warping functions
%   setup                         - setup environment
%   setup_paths                   - setup paths
%   ShiftF                        - Shifts curve starting point by tau
%   simul_align                   - Align two functions by simultaneous reparameterization
%   simul_gam                     - Align two function using simultaneous alignments.
%   smooth_data                   - Smooth Functions
%   SqrtMean                      - SRVF transform of warping functions calculates mean
%   SqrtMeanInverse               - SRVF transform of warping functions calculates mean
%   SqrtMedian                    - SRVF transform of warping functions calculates median
%   srvf_to_f                     - Convert SRSF to f
%   Translation_Boxplot           - Translation Boxplot
%   v_to_gam                      - tangent space to warping function
%   warp_curve_gamma              - Warp curve by gamma
%   warp_f_gamma                  - Warp Function by gamma
%   warp_q_gamma                  - WARP_F_GAMMA Warp SRVF by gamma

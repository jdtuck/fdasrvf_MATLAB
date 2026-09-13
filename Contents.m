% FDASRVF
%
% Files
%   align_fPCA                       - Group-wise function alignment and PCA Extractions
%   ampbox                           - Produce amplitude box plots
%   bootTB                           - Bootstrapped functional tolerance bounds
%   calculateCentroid                - Calculate centroid of curve
%   create_basismatrix               - Create B-spline basis matrix
%   curve_to_q                       - Convert curve to Square-Root Velocity Function
%   cycl_procrustes                  - Align two matrices by a cyclic row shift and an orthogonal transform
%   cycl_procrustes_fft              - Fast cyclic Procrustes using FFT
%   elastic_depth                    - Compute elastic depth
%   elastic_distance                 - Calculates the two elastic distances between two functions
%   elastic_distance_curve           - Calculates the elastic distance between two curves
%   elastic_logistic                 - A class to provide SRVF logistic regression
%   elastic_lpcr_regression          - A class to provide SRVF logistic principal component regression
%   elastic_mlogistic                - A class to provide SRVF multinomial logistic regression
%   elastic_mlpcr_regression         - A class to provide SRVF multinomial logistic principal component regression
%   elastic_pcr_regression           - A class to provide SRVF principal component regression
%   elastic_regression               - A class to provide SRVF regression
%   exp_map                          - Exponential Map
%   f_to_srvf                        - Convert function to Square-Root Velocity Function
%   fdacurve                         - A class to provide registration of curves in R^n using SRVF
%   fdahpca                          - A class to provide horizontal fPCA
%   fdahpns                          - A class to provide horizontal fPNS
%   fdajpca                          - A class to provide joint fPCA
%   fdakma                           - A class to provide a kmeans clustering and alignment
%   fdavpca                          - A class to provide vertical fPCA
%   fdawarp                          - A class to provide alignment of functional data using SRVF
%   gam_to_h                         - Convert warping function to h space
%   gam_to_psi                       - Convert warping function to hilbert sphere
%   gam_to_v                         - Convert warping function to tangent space on hilbert sphere
%   geodesic_sphere_Full             - Calculates geodesic on sphere
%   h_to_gam                         - Convert h to warping function
%   inv_exp_map                      - Inverse Exponential Map
%   invertGamma                      - Invert Warping Function
%   joint_gauss_model                - Gaussian generative model from joint fPCA
%   L2norm                           - L2 Functional Norm
%   optimum_reparam                  - Calculates Warping for two SRVFs
%   optimum_reparam_curve            - Calculates warping for two curve SRVFs
%   outlier_detection                - Outlier Detection
%   pairwise_align                   - Align two functions
%   pairwise_align_bayes             - Align two functions using Bayesian method
%   pairwise_align_bayes_infHMCbasis - Align two trajectories using hierarchical Bayesian HMC
%   pairwise_align_curves            - Registers two curves
%   Path_Plot                        - Plot geodesic path
%   pcaTB                            - Functional principal component tolerance bound generation
%   phbox                            - Construct phase box plots
%   plot_curve                       - Plot a curve in R^2
%   PNS_warping                      - Principal nested spheres of warping functions
%   psi_to_gam                       - Hilbert sphere to warping function
%   q_to_curve                       - Convert SRVF to curve
%   randomGamma                      - Generative model for warping functions
%   ReSampleCurve                    - Resample curve to have N points
%   rgam                             - Generate random warping functions
%   setup_paths                      - Setup paths
%   ShiftF                           - Shifts curve starting point by tau
%   simul_align                      - Align two functions by simultaneous reparameterization
%   simul_gam                        - Align two function using simultaneous alignments
%   smooth_data                      - Smooth Functions
%   SqrtMean                         - SRVF transform of warping functions calculates mean
%   SqrtMeanInverse                  - SRVF transform of warping functions calculates inverse mean
%   SqrtMedian                       - SRVF transform of warping functions calculates median
%   srvf_to_f                        - Convert SRSF to f
%   Translation_Boxplot              - Translation Boxplot
%   v_to_gam                         - Tangent space to warping function
%   warp_curve_gamma                 - Warp curve by gamma
%   warp_f_gamma                     - Warp Function by gamma
%   warp_q_gamma                     - Warp SRVF by gamma
%   warp_srvf_gamma                  - Warp SRVF by gamma (curves)

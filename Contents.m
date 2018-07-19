% FDASRVF
%
% Files
%   align_fPCA                    - Group-wise function alignment and PCA Extractions
%   ampbox                        - ampbox Produce amplitude box plots
%   Basis_Normal_A                - Find Normal Basis
%   bootstrap_tol                 - Boostrapped tolerance bounds
%   calculateCentroid             - Calculate centroid of curve
%   create_basismatrix            - Create Basis Matrix
%   cumsimps                      - Cumulative Simpson's numerical integration.
%   cumtrapzmid                   - Cummulative Trapezodial Integration from mid point
%   curve_to_q                    - Convert curve to Square-Root Velocity Function
%   elastic_distance              - Calculates the two elastic distances between two
%   elastic_logistic              - elastic_logistic A class to provide SRVF logistic regression
%   elastic_mlogistic             - elastic_mlogistic A class to provide wSRVF multionomial logistic
%   elastic_regression            - elastic_regression A class to provide SRVF regression
%   exp_map                       - Exponential Map
%   f_to_srvf                     - Convert function to Square-Root Velocity Function
%   fdacurve                      - fdacurve A class to provide registration of curves in R^n using SRVF
%   fdahpca                       - fdahpca A class to provide horizontal fPCA
%   fdakma                        - fdakma A class to provide a kmeans clustering and alignment
%   fdavpca                       - fdavpca A class to provide vertical fPCA
%   fdawarp                       - fdawarp A class to provide alignment of functional data using SRVF
%   Find_Best_Rotation            - Find best rotation between two curves
%   Find_Rotation_and_Seed_unique - find roation and seed of two curves
%   geodesic_sphere_Full          - Calculates geodesic on sphere
%   Gram_Schmidt                  - orthonormalize basis
%   inner_product                 - Functional inner product
%   InnerProd_Q                   - Calculate inner product of two SRVFs
%   inv_exp_map                   - Inverse Exponential Map
%   invertGamma                   - Invert Warping Function
%   L2norm                        - L2 Functional Norm
%   optimum_reparam               - Calculates Warping for two SRVFs
%   outlier_detection             - Outlier Detection
%   pairwise_align                - Align two functions
%   pairwise_align_bayes          - PAIRWISE_ALIGN Align two functions using Bayesian method
%   pairwise_align_curves         - registers to curves
%   Parallel_Transport_C          - PARALLEL_TRANSPORT Parallel transport vector
%   Path_Plot                     - Plot geodesic path
%   phbox                         - phbox Construct phase box plots
%   progress                      - Wrapper class to provide an iterator object for loop creation
%   ProgressBar                   - A class to provide a convenient and useful progress bar
%   Project_Tangent               - Project tangent
%   ProjectC                      - Project closed curve
%   q_to_curve                    - Convert SRVF to curve
%   ReSampleCurve                 - Resample curve to have N points
%   rgam                          - Generate random warping functions
%   setup                         - setup environment
%   setup_paths                   - setup paths
%   ShiftF                        - Shifts curve starting point by tau
%   simps                         - Simpson's numerical integration.
%   simul_align                   - Align two functions by simultaneous reparameterization
%   simul_gam                     - Align two function using simultaneous alignments.
%   smooth_data                   - Smooth Functions
%   SqrtMean                      - SRVF transform of warping functions calculates mean
%   SqrtMeanInverse               - SRVF transform of warping functions calculates mean
%   SqrtMedian                    - SRVF transform of warping functions calculates median
%   srvf_to_f                     - Convert SRSF to f
%   Translation_Boxplot           - Translation Boxplot
%   updateParallel                - Update function when ProgressBar is used in parallel setup
%   varycolor                     - Produces colors with maximum variation on plots with multiple
%   warp_curve_gamma              - Warp curve by gamma
%   warp_f_gamma                  - Warp Function by gamma
%   warp_q_gamma                  - WARP_F_GAMMA Warp SRVF by gamma

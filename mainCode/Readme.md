Code for implementing mtclust

The main function for implementing mtclust is m_t_clust. Note that the
code is somehow commented. In your use it is recommended that several
starting solutions are compared. Occasional error messages might be
obtained if the starting solution is particularly bad. 

The m_t_clust function takes in input:  

- data: a data matrix or a dataframe with numerical entries with observations
  on the rows and the variables on the columns

- alpha.fixed: the desired trimming level

- K: the number of clusters

- modelNames: the desired parametrization of the covariance
  matrix. Refer to mclust R package vignettes 
for the different options available

- restr.factor: the restriction factor to be imposed on the eigenvalue ratio 
  (needed only for VVI, EVE, EVV and VVV models)

- alpha.init: the model is initialized with the tclust model with alpha.init as initial trimming level

- restr.fact.init: the initial restriction factor for tclust initialization. 
  Defaults to 50. By combining restr.fact.init=1 and alpha.init=0
  the standard k-means initialization is obtained

- n.start.init: number of random starts in the initialization

- initial.values: a list with the initial values that can be used to
   by-pass the random or deterministic initializations above. 
  The list of values must contain the following elements:

  1) $cov 

  an array whose dimensions are p,p,K where p is the number of variables and
  K is the number of desired clusters. Each slide of the array
  contains the covariance
  matrix of each cluster. 

 2) $sigmasq
 A vector whose length must be equal to K containing the p-th  oot of the deteriminant of each
 of the scatter matrix

 3) $shape
 a matrix with p rows and K columns. 
Each column has to contain the normalized eigenvalues 
 of each covariance matrix

 4) $orientation
 an array with dimensions p,p,K. Each slide of the array contains a
 matrix whose columns are given by the eigenvetors of each covariance matrix

 5) $cholsigma
 an array containing the Cholesky decomposition of each scatter
 matrix. This can be easily obtained though the function chol()

 6) $weights
 a vector of length K with cluster proportions

 7) $mean 
 a matrix whose columns are the vector mean of each cluster

 8) $cluster
 a vector whose length must be equal to number of the observations,
 with the initial clustering. Trimmed values 
 must be labelled as 0

The output is as: 

 1) $data: the dataset used for the clustering model

 2) $classification: The estimated classification labels. Outlying observations are
    labeled as 0

 3) $init the partition provided for the initialization of the algorithm. If no intialization
    has been provided by the user the initialization with its default values is returned

 4) $loglikelihood A vector containing the values of the loglikelihood function at each iteration of the CEM algorithm

 5) $contribution a vector containing the contribution of the likelihood function of each observation,     including outlying observations

 6) $parameters: a list containing the following elements:
     $pro: the estimated clusters proportions
     $mean: A matrix continaining in each column the vector mean of each cluster
     $variance: A list containing all the variance components. Depending on 
     the parametrization it will contain different elements. 

 7) $war: a logical value, indicating whether the ER has been actually
 used to restrict the covariance estimates during the CEM algorithm iterations

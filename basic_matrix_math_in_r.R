# Fish 507: Applied Time-Series Analysis in Fisheries and Environmental Sciences
# Basic matrix math in R

## 1.1 Creating matrices in R
## Create a 3 × 4 matrix, meaning 3 row and 4 columns, that is all 1s:
matrix(data = 1, nrow = 3, ncol = 4)
## OUTPUT:
#       [,1] [,2] [,3] [,4]
# [1,]    1    1    1    1
# [2,]    1    1    1    1
# [3,]    1    1    1    1

## Create a 3×4 matrix filled in with the numbers 1 to 12 by column (default) and by row:
matrix(1:12,3,4)
## OUTPUT:
#       [,1] [,2] [,3] [,4]
# [1,]    1    4    7   10
# [2,]    2    5    8   11
# [3,]    3    6    9   12
matrix(1:12,3,4, byrow = TRUE)
## OUTPUT:
#       [,1] [,2] [,3] [,4]
# [1,]    1    2    3    4
# [2,]    5    6    7    8
# [3,]    9   10   11   12

## Create a matrix with one column:
matrix(1:6, ncol = 1)
## OUTPUT:
#      [,1]
# [1,]    1
# [2,]    2
# [3,]    3
# [4,]    4
# [5,]    5
# [6,]    6

## Create a matrix with one row:
matrix(1:6, nrow = 1)
## OUTPUT:
#       [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]    1    2    3    4    5    6

## Check the dimensions of a matrix:
foo = matrix(1:6,3,3)
dim(foo)
## OUTPUT: 
# [1] 3 3

## Get the number of rows in a matrix:
nrow(foo)
## OUTPUT:
# [1] 3
## ALTERNATE METHOD:
dim(foo)[1]
## OUTPUT:
# [1] 3

## Create a 3D matrix (called array):  
A=array(1:6, dim=c(2,3,2)) # translation: make a matrix(1:6,2,3)  2 times and store in an array.
A
## OUTPUT:
# , , 1
# 
#       [,1] [,2] [,3]
# [1,]    1    3    5
# [2,]    2    4    6
# 
# , , 2
# 
#       [,1] [,2] [,3]
# [1,]    1    3    5
# [2,]    2    4    6
## Check dimensions of array:
dim(A)
## OUTPUT:
# [1] 2 3 2

## Check if an object is a matrix. A dataframe is not a matrix. A vector is not a matrix.
A=matrix(1:4, 1, 4)
A
## OUTPUT:
#       [,1] [,2] [,3] [,4]
# [1,]    1    2    3    4
class(A)
## OUTPUT:
# [1] "matrix"
B=data.frame(A)
B
## OUTPUT:
#   X1 X2 X3 X4 
# 1  1  2  3  4
class(B)
## OUTPUT:
# [1] "data.frame"
C=1:4 
C
## OUTPUT:
# [1] 1 2 3 4
class(C)
## OUTPUT:
# [1] "integer"

## 1.2 Matrix multiplication, addition and transpose
## In R, you use the %*% operation to do matrix multiplication. When you do matrix multiplication,  ##
## the columns of the matrix on the left must equal the rows of the matrix on the right. The result ##
## is a matrix that has the number of rows of the matrix on the left and number of columns of the   ##
## matrix on the right.
# (n×m)(m× p) = (n× p)
A=matrix(1:6, 2, 3) #2 rows, 3 columns
B=matrix(1:6, 3, 2) #3 rows, 2 columns
A%*%B #this works
## OUTPUT:
#       [,1] [,2]
# [1,]   22   49
# [2,]   28   64
B%*%A #this works
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    9   19   29
# [2,]   12   26   40
# [3,]   15   33   51
try(B%*%B)
## OUTPUT:
# Error in B %*% B : non-conformable arguments

## To add two matrices use +. The matrices have to have the same dimensions.
A+A #works
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    2    6   10
# [2,]    4    8   12
A+t(B) #works
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    2    5    8
# [2,]    6    9   12
try(A+B) #does not work since A has 2 rows and B has 3
## OUTPUT:
# Error in A + B : non-conformable arrays

## The transpose of a matrix is denoted A⊤ or A′. To transpose a matrix in R, you use t().
A=matrix(1:6, 2, 3) #2 rows, 3 columns
t(A)
## OUTPUT:
#       [,1] [,2]
# [1,]    1    2
# [2,]    3    4
# [3,]    5    6

## 1.3 Subsetting a matrix
## To subset a matrix, we use [ ]:
A=matrix(1:9, 3, 3) #3 rows, 3 columns 
## Get the first and second rows of A: 
## NOTE: this result will be a 2x3 matrix
A[1:2,]
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    1    3    5
# [2,]    2    4    6

## Get the top 2 rows and left 2 columns:
A[1:2,1:2]
## OUTPUT:
#       [,1] [,2]
# [1,]    1    4
# [2,]    2    5

## What does this do?
A[c(1,3),c(1,3)] # get the first and third row entries, first and third column entries
## OUTPUT:
#       [,1] [,2]
# [1,]    1    7
# [2,]    3    9

## This?
A[c(1,2,1),c(2,3)] # get the first, second, first again row entries from the second and third columns
## OUTPUT:
#       [,1] [,2]
# [1,]    4    7
# [2,]    5    8
# [3,]    4    7

## Get the element of a matrix in row 1 and the last column:
A=matrix(1:9, 3, 3)
A[1,ncol(A)]
## OUTPUT:
# [1] 7
## ALTERNATE METHOD:
A[1,dim(A)[2]]
## OUTPUT:
# [1] 7

## IMPORTANT NOTE: R will create vectors from subsetting matrices! Happens with subsets with 1 row or 1 column
## SOLUTION: use drop=FALSE
## Example:
A=matrix(1:9, 3, 3)
## Get the first row:
B=A[1,,drop=FALSE]

## 1.4 Replacing elements in a matrix
## Replace 1 element:
A=matrix(1, 3, 3)
A[1,1]=2
A
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    2    1    1
# [2,]    1    1    1
# [3,]    1    1    1

## Replace a row with all 2s:
A=matrix(1, 3, 3)
A[1,]=2
A
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    2    2    2
# [2,]    1    1    1
# [3,]    1    1    1

## Replace a row with numbers 1:3:
A[1,]=1:3 
A
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    1    2    3
# [2,]    1    1    1
# [3,]    1    1    1

## Replace group of elements. This often does not work as one expects so be sure look at  ##
## your matrix after trying something like this. Here I want to replace elements (1,3) and##
## (3,1) with 2, but it didn’t work as I wanted.
A=matrix(1, 3, 3)
A[c(1,3),c(3,1)]=2
A
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    2    1    2
# [2,]    1    1    1
# [3,]    2    1    2

## How do I replace elements (1,1) and (3,3) with 2 then? It’s tedious. If you have a lot##
## of elements to replace, you might want to use a for loop.
A=matrix(1, 3, 3)
A[1,3]=2
A[3,1]=2
A
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    1    1    2
# [2,]    1    1    1
# [3,]    2    1    1

## 1.5 Diagonal matrices and identity matrices
## A diagonal matrix is one that is square, meaning number of rows equals num- ber of   ##
## columns, and it has 0s on the off-diagonal and non-zeros on the diagonal. In R, you  ##
## form a diagonal matrix with the diag() function:
diag(1,3) #put 1 on diagonal of 3x3 matrix
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    1    0    0
# [2,]    0    1    0
# [3,]    0    0    1

diag(1:4) #put 1 to 4 on diagonal of 4x4 matrix
## OUTPUT:
#       [,1] [,2] [,3] [,4]
# [1,]    1    0    0    0
# [2,]    0    2    0    0
# [3,]    0    0    3    0
# [4,]    0    0    0    4

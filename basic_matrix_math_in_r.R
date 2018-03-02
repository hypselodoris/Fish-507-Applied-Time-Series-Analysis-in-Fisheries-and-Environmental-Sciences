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

## The diag() function can also be used to replace elements on the diagonal of a matrix:
A=matrix(3, 3, 3)
diag(A)=1
A
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    1    3    3
# [2,]    3    1    3
# [3,]    3    3    1

A=matrix(3, 3, 3)
diag(A)=1:3
A
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    1    3    3
# [2,]    3    2    3
# [3,]    3    3    3

A=matrix(3, 3, 4)
diag(A[1:3,2:4])=1
A
## OUTPUT:
#       [,1] [,2] [,3] [,4]
# [1,]    3    1    3    3
# [2,]    3    3    1    3
# [3,]    3    3    3    1

## The diag function is also used to get the diagonal of a matrix:
A=matrix(1:9, 3, 3)
diag(A)
## OUTPUT:
# [1] 1 5 9

## The identity matrix is a special kind of diagonal matrix with 1s on the diagonal.  ##
## It is denoted I. I3 would mean a 3×3 diagonal matrix. A identity matrix has the    ##
## property that AI=A and IA=A so it is like a 1.
A=matrix(1:9, 3, 3)
I=diag(3) #shortcut for 3x3 identity matrix
A%*%I
## OUTPUT:
#       [,1] [,2] [,3]
# [1,]    1    4    7
# [2,]    2    5    8
# [3,]    3    6    9

## 1.6 Taking the inverse of a square matrix                                          ##
## The inverse of a matrix is denoted A−1. You can think of the inverse of a matrix   ##
## like 1/a. 1/a×a = 1. A−1A = AA−1 = I. The inverse of a matrix does not always      ##
## exist; for one it has to be square. We’ll be using inverses for variance-covariance##
## matrices and by definition (of a variance-covariance matrix), the inverse of those ##
## exist. In R, there are a couple way common ways to take the inverse of a variance- ##
## covariance matrix (or something with the same properties). solve is the most common##
## probably:
A=diag(3,3)+matrix(1,3,3)
invA=solve(A)
invA%*%A
## OUTPUT:
#               [,1]          [,2] [,3]
# [1,] 1.000000e+00 -6.938894e-18    0
# [2,] 2.081668e-17  1.000000e+00    0
# [3,] 0.000000e+00  0.000000e+00    1

A%*%invA
## OUTPUT:
#               [,1]          [,2] [,3]
# [1,] 1.000000e+00 -6.938894e-18    0
# [2,] 2.081668e-17  1.000000e+00    0
# [3,] 0.000000e+00  0.000000e+00    1

## Another option is to use chol2inv which uses a Cholesky decomposition:
## NOTE: The Cholesky decomposition is a handy way to keep your variance-covariance ##
## matrices valid when doing a parameter search. Don’t search over the raw variance-##
## covariance matrix. Search over a matrix where the lower triangle is 0, that is   ##
## what a Cholesky decomposition looks like. Let’s call it B. Your variance-        ##
## covariance matrix is t(B)%*%B.                                                   ##
## NOTE: For the purpose of this course, solve is fine.                             ##
A=diag(3,3)+matrix(1,3,3)
invA=chol2inv(chol(A))
invA%*%A
## OUTPUT:
#               [,1]         [,2]          [,3]
# [1,]  1.000000e+00 6.938894e-17  0.000000e+00
# [2,]  2.081668e-17 1.000000e+00 -2.775558e-17
# [3,] -5.551115e-17 0.000000e+00  1.000000e+00
A%*%invA
## OUTPUT:
#               [,1]          [,2]          [,3]
# [1,] 1.000000e+00  2.081668e-17 -5.551115e-17
# [2,] 6.938894e-17  1.000000e+00  0.000000e+00
# [3,] 0.000000e+00 -2.775558e-17  1.000000e+00

## PROBLEMS ##
## 1.1 Build a 4×3 matrix with the numbers 1 through 4 in each row. ##
A=matrix(1:4,nrow = 4, ncol = 3)

## 1.2 Extract the elements in the 1st and 2nd rows and 1st and 2nd columns ##
## (you’ll have a 2×2 matrix). Show the R code that will do this.           ##
B=A[1:2, 1:2]
B
## OUTPUT:
#       [,1] [,2]
# [1,]    1    1
# [2,]    2    2

## 1.3 Build a 4×3 matrix with the numbers 1 through 12 by row (meaning the ##
## first row will have the numbers 1 through 4 in it).                      ##
C=matrix(1:12, 4, 3, byrow = TRUE)

## 1.4 Extract the 3rd row of the above. Show R code to do this where you   ##
## end up with a vector and how to do this where you end up with a 1×3      ##
## matrix.
# Vector result:
D=C[3,]
# Matrix result:
D=C[3,,drop=FALSE]

## 1.5 Build a 4×3 matrix that is all 1s except a 2 in the (2,3) element    ##
## (2nd row, 3rd column).
A = matrix(1,4,3)
A[2,3]=2

## 1.6 Take the transpose of the above.                                     ##
t(A)

## 1.7 Build a 4 × 4 diagonal matrix with 1 through 4 on the diagonal.      ##
A = diag(1:4)

## 1.8 Build a 5×5 identity matrix.
I=diag(5)

## 1.9 Replace the diagonal in the above matrix with 2 (the number 2).      ##
diag(I)=2

## 1.10 Build a matrix with 2 on the diagonal and 1s on the offdiagonals.   ##
A = matrix(1,2,2)
A[1,1]=2
A[2,2]=2
# More elegant solution:
A = diag(1,4)+1

## 1.11 Take the inverse of the above.
invA = solve(A)

## 1.12 Build a 3×3 matrix with the first 9 letters of the alphabet. First  ##
## column should be “a”, “b”, “c”. letters[1:9] gives you these letters.
A = matrix(letters[1:9],3,3)

## 1.13 Replace the diagonal of this matrix with the word “cat”.            ##
diag(A) <- c("c", "a", "t")

## 1.14 Build a 4×3 matrix with all 1s. Multiply by a 3×4 matrix with all 2s##
A = matrix(1,4,3)
B = matrix(2,3,4)
A%*%B

## 1.15 If A is a 4×3 matrix, is AA possible? Is AA> possible? Show how to ##
## write AA> in R.
A = matrix(2,4,3)
A%*%A # This throws error: Error in A %*% A : non-conformable arguments
A%*%t(A) # This works

## 1.16 In the equation, AB = C, let A =                                  ##
#       [,1] [,2] [,3]
# [1,]    1    4    7
# [2,]    2    5    8
# [3,]    3    6    9
## Build a B matrix with only 1s and 0s such that the values on the       ##
## diagonal of C are 1, 8, 6 (in that order).Show your R code for A, B    ##
## and AB.
A = matrix(1:9,3,3)
B = matrix(0,3,3)
B[1,1] = 1
B[3,2] = 1
B[2,3] = 1
C = A%*%B
diag(C)

## 1.17 Same A matrix as above and same equation AB = C. Build a 3 × 3 B  ##
## matrix such that C = 2A. So C =
# h
# 2 8 14
# 4 10 16
# 6 12 18
# i
# . Hint, B is diagonal.
A = matrix(1:9,3,3)
C = 2*A
B = diag(2,3)
A%*%B

## 1.18 Same A and AB = C equation. Build a B matrix to compute the row   ##
## sums of A. So the first ‘row sum’ would be 1+4+7, the sum of all       ##
## elements in row 1 of A. C will be 
# h
# 12
# 15
# 18
# i
# , the row sums of A. Hint, B is a column matrix
# (1 column).
A = matrix(1:9,3,3)
B = matrix(1,3,1)
A%*%B

## 1.19 Same A matrix as above but now equation BA = C. Build a B       ##
## matrix to compute the column sums of A. So the first ‘column sum’    ##
## would be 1+2+3. C will be a 1×3 matrix.
A = matrix(1:9,3,3)
B = matrix(1,1,3)
B%*%A

## 1.20 Let AB = C equation but A =
# h
# 2 1 1
# 1 2 1
# 1 1 2 i
# (so A=diag(3)+1). Build a B matrix
# such that C =
#   h
# 3
# 3
# 3
# i
# . Hint, you need to use the inverse of A.
A = diag(3)+1
C=matrix(3,3,1)
#AB=C
#B=inv(A)%*%C
B=solve(A)%*%C
B

(* :Title: Matrix.m *)

(* :Author: Mauricio C. de Oliveira *)

(* :Context: Matrix` *)

(* :Sumary: *)

(* :Alias: *)

(* :Warnings: *)

(* :History:
   :3/13/13: Deprecated LinearAlgebra`MatrixManipulation`. (mauricio)
   :06/01/2004: First Version
*)

BeginPackage[ 
     "Matrix`",
     "NCLinearAlgebra`"
];

Clear[Matrix]

Matrix::usage = "
Matrix[m] wraps a two-dimensional list under a Matrix head. \ 
Try to substitute 'a + b.c' with matrices by applying a rule for one \
matrix at a time and see what happens... \ 
This is all because Plus is Listable and Mathematica makes no \
distinction between two-dimensional lists and matrices.";


Clear[MatrixToList]

MatrixToList::usage = "MatrixToList[m] converts Matrix m into a two dimensional list.";


Clear[MatrixToVector]

MatrixToList::usage = "MatrixToVector[m] converts Matrix m into a one dimensional list (column oriented).";


Clear[SymmetricMatrixToVector]

MatrixToList::usage = "SymmetricMatrixToVector[m] converts symmetric Matrix m into a one dimensional list (lower diagonal, column oriented).";


Clear[MatrixFlatten]

MatrixFlatten::usage = "MatrixFlatten[m] attempts to flatten first level entries in a Matrix. If it fails it returns the original Matrix";

Clear[NumericMatrixQ]

NumericMatrixQ::usage = "NumericMatrixQ[m] returs True if all entries \
of Matrix return True to NumericQ";

Clear[ZeroQ]

ZeroQ::usage = "ZeroQ[m] returs True if all entries \
of Matrix are zero";

Clear[VectorToMatrix]

VectorToMatrix::usage = "VectorToMatrix[m] converts a Mathematica vector (unidimensional list) into a column Matrix.";

Clear[MatrixDimensions]

MatrixDimensions::usage = "MatrixDimensions[m] provides Dimensions that accounts for Matrix entries.";

Clear[MatrixProductDimensions]
Clear[MatrixEntryDimensions]


Begin["`Private`"]; 

  (* Messages *)

  Matrix::NotMatrix = "This is not a Matrix.";
  Matrix::WrongDimensions="Incompatible dimensions."


  (* Attributes *)

  SetAttributes[Matrix, NumericFunction];


  (* Constructors *)

  Matrix[x_?MatrixQ] := Apply[Matrix, x];
  Matrix[x___] := Matrix[] /; !MatrixQ[{x}];
  Matrix[] := (Message[Matrix::NotMatrix]; $Failed);

  (* Conversions *)
  
  MatrixToList[a_Matrix, Infinity] := ReplaceRepeated[a, x_Matrix->MatrixToList[x]];
  MatrixToList[a_Matrix] := Apply[List,a];

  (* Part *)

  Clear[MatrixPartAux];
  MatrixPartAux[x_?MatrixQ] := Matrix[x];
  MatrixPartAux[x_] := x;
  
  Part[m_Matrix, i_, j_] ^:= MatrixPartAux[Part[MatrixToList[m], i, j]];

  (* Dimensions *)

  Dimensions[m_Matrix] ^:= Dimensions[MatrixToList[m]];
  
  (* ArrayFlatten *)

  ArrayFlatten[a_Matrix] ^:= Matrix[ArrayFlatten[MatrixToList[a,Infinity]]];
  
  (* Numeric tests *)

  NumberQ[a_Matrix] ^:= MatrixQ[MatrixToList[a], NumberQ];

  NumericMatrixQ[a_Matrix] := MatrixQ[MatrixToList[a], NumericQ];
  NumericMatrixQ[a_] := False;

  ZeroQ[a_Matrix] ^:= Apply[And, Flatten[Map[PossibleZeroQ, MatrixToList[a], {2}]]];

  (* Ouput *)

  MatrixForm[a_Matrix] ^:= MatrixForm[MatrixToList[a]]

  Format[x_Matrix, StandardForm] ^:= MatrixForm[x]

  
  (* Arithmetic Operations *)

  Plus[a_Matrix, b_Matrix] ^:= Matrix[Plus[MatrixToList[a], MatrixToList[b]]];

  Times[a_Matrix, b_Matrix] ^:= Matrix[Times[MatrixToList[a], MatrixToList[b]]];
  Times[a_, b_Matrix] ^:= Matrix[a MatrixToList[b]]

  Power[a_Matrix, b_] ^:= Matrix[Power[MatrixToList[a], b]];


  (* Matrix Algebra *)

  Dot[a_Matrix, b_Matrix] ^:= Matrix[Dot[MatrixToList[a], MatrixToList[b]]];

  Transpose[a_Matrix] ^:= Matrix[Transpose[MatrixToList[a]]];

  Conjugate[a_Matrix] ^:= Matrix[Conjugate[MatrixToList[a]]];

  Inverse[a_Matrix] ^:= Matrix[Inverse[MatrixToList[a]]];

  Det[a_Matrix] ^:= Det[MatrixToList[a]];
  Det[a_Matrix, option_] ^:= Det[MatrixToList[a], option];

  Tr[a_Matrix] ^:= Tr[MatrixToList[a]];
  Tr[a_Matrix, f_] ^:= Tr[MatrixToList[a], f];
  Tr[a_Matrix, f_, n_] ^:= Tr[MatrixToList[a], f, n];


  (* Linear Algebra *)

  Eigenvalues[a_Matrix] ^:= Eigenvalues[MatrixToList[a]];
  Eigenvalues[a_Matrix, b_Matrix] ^:= Eigenvalues[{MatrixToList[a], MatrixToList[b]}];
  Eigenvalues[a_Matrix, k_] ^:= Eigenvalues[MatrixToList[a], k];

  Eigenvectors[a_Matrix] ^:= Eigenvectors[MatrixToList[a]];
  Eigenvectors[a_Matrix, b_Matrix] ^:= Eigenvectors[{MatrixToList[a], MatrixToList[b]}];
  Eigenvectors[a_Matrix, k_] ^:= Eigenvectors[MatrixToList[a], k];

  Eigensystem[a_Matrix] ^:= Eigensystem[MatrixToList[a]];
  Eigensystem[a_Matrix, b_Matrix] ^:= Eigensystem[{MatrixToList[a], MatrixToList[b]}];
  Eigensystem[a_Matrix, k_] ^:= Eigensystem[MatrixToList[a], k];

  CharacteristicPolynomial[a_Matrix, x_] ^:= CharacteristicPolynomial[MatrixToList[a], x];

  LinearSolve[a_Matrix, b_Matrix] ^:= Matrix[LinearSolve[MatrixToList[a], MatrixToList[b]]];

  NullSpace[a_Matrix] ^:= Matrix[NullSpace[MatrixToList[a]]];

  MatrixRank[a_Matrix] ^:= MatrixRank[MatrixToList[a]];

  RowReduce[a_Matrix] ^:= Matrix[RowReduce[MatrixToList[a]]];

  Minors[a_Matrix] ^:= Matrix[Minors[MatrixToList[a]]];
  Minors[a_Matrix, k_] ^:= Matrix[Minors[MatrixToList[a], k]];

  MatrixPower[a_Matrix, n_] ^:= Matrix[MatrixPower[MatrixToList[a], n]];

  MatrixExp[a_Matrix] ^:= Matrix[MatrixExp[MatrixToList[a]]];

  (* Outer[f_, a_Matrix] ^:= MatrixExp[MatrixToList[a]]; *)

  Norm[a_Matrix] ^:= Norm[MatrixToList[a]];
  Norm[a_Matrix, p_] ^:= Norm[MatrixToList[a], p];

  SingularValueList[a_Matrix] ^:= SingularValueList[MatrixToList[a]];
  SingularValueList[a_Matrix, b_Matrix] ^:= SingularValueList[{MatrixToList[a], MatrixToList[b]}];
  SingularValueList[a_Matrix, k_] ^:= SingularValueList[MatrixToList[a], k];

  SingularValueDecomposition[a_Matrix] ^:= Map[Matrix, SingularValueDecomposition[MatrixToList[a]]];
  SingularValueDecomposition[a_Matrix, b_Matrix] ^:= Map[Matrix, SingularValueDecomposition[{MatrixToList[a], MatrixToList[b]}] ];
  SingularValueDecomposition[a_Matrix, k_] ^:= Map[Matrix, SingularValueDecomposition[MatrixToList[a], k] ];

  PseudoInverse[a_Matrix] ^:= PseudoInverse[MatrixToList[a]];

  QRDecomposition[a_Matrix] ^:= Map[Matrix, QRDecomposition[MatrixToList[a]]];

  LUDecomposition[a_Matrix] ^:= MapAt[Matrix, LUDecomposition[MatrixToList[a]], {1}];

  CholeskyDecomposition[a_Matrix] ^:= Matrix[CholeskyDecomposition[MatrixToList[a]]];

  SchurDecomposition[a_Matrix] ^:= Map[Matrix, SchurDecomposition[MatrixToList[a]]];
  SchurDecomposition[a_Matrix, b_Matrix] ^:= Map[Matrix, SchurDecomposition[{MatrixToList[a], MatrixToList[b]}]];

  JordanDecomposition[a_Matrix] ^:= Map[Matrix, JordanDecomposition[MatrixToList[a]]];
  
  



End[];

EndPackage[];

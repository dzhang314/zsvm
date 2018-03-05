(* ::Package:: *)

packedMatrix[a_String, n_Integer?Positive] := Block[
	{x = ConstantArray[0, {n, n}], k = 0},
	Do[x[[i, j]] = a[k++], {i, n}, {j, i}];
	x + Transpose@LowerTriangularize[x, -1]];

negativeTermQ[term_] := (Head[term] === Times) &&
	AnyTrue[term, NumericQ[#] && Negative[#]&];
coercePositiveTerm[term_] := If[
	negativeTermQ[term], -term, term, term];

toCCode[n_Integer] := ToString[n];
toCCode[x_String[i_Integer]] := x <> "[" <> toCCode[i] <> "]";
toCCode[Power[x_, 2]] := With[{s = toCCode[x]}, s <> " * " <> s];

toCCode[Times[x_, xs__]] := StringJoin@Riffle[
	toCCode /@ {x, xs}, " * "];
toCCode[Plus[x_, xs__]] := Module[
	{signPattern = negativeTermQ /@ {x, xs}},
	StringJoin@Riffle[Prepend[
		<|True -> "\n- ", False -> "\n+ "|> /@ Rest[signPattern],
		<|True -> "-", False -> ""|>@First[signPattern]],
		toCCode @* coercePositiveTerm /@ {x, xs}]];

determinantCCode[n_Integer?Positive] :=
	"const double det =\n" <>
		toCCode@Expand@Det@packedMatrix["x", n] <> ";\n";

expandedAdjugateEntries[n_Integer?Positive] := With[
	{x = packedMatrix["x", n]},
	toCCode /@ Expand[Join @@ MapThread[Take, {
		LowerTriangularize[Det[x] * Inverse[x]],
		Range[n]}]]];

adjugateCCode[n_Integer?Positive] :=
	StringJoin@MapThread[StringJoin, {
		"y[" <> ToString[#] <> "] = ("& /@
			Range[0, n * (n + 1) / 2 - 1],
		expandedAdjugateEntries[n],
		ConstantArray[") / det;\n", n * (n + 1) / 2]}];

determinantInverseFunctionCCode[n_Integer?Positive]:=StringJoin[
"double packed_determinant_inverse_",ToString[n],
"(\nconst double *__restrict__ x, double *__restrict__ y) {\n",
determinantCCode[n],adjugateCCode[n],"return det;\n}\n"];

kineticTraceFunctionCCode[n_Integer?Positive]:=
"double packed_kinetic_trace_"<>ToString[n]<>"(
const double *__restrict__ a, const double *__restrict__ b,
const double *__restrict__ c, const double *__restrict__ m) {
return "<>toCCode@Expand@Tr@Dot[
packedMatrix["a",n],packedMatrix["c",n],
packedMatrix["b",n],DiagonalMatrix["m"/@Range[0,n-1]]
]<>";\n}\n";

quadraticFormFunctionCCode[n_Integer?Positive]:=
"double packed_quadratic_form_"<>ToString[n]<>"(
const double *__restrict__ x, const double *__restrict__ v) {
return "<>toCCode@Expand@Dot["v"/@Range[0,n-1],
packedMatrix["x",n],"v"/@Range[0,n-1]]<>";\n}\n";

permutationConjugateFunctionCCode[n_Integer?Positive]:=
With[{p=ArrayReshape["p"/@Range[0,n^2-1],{n,n}]},
StringJoin["void packed_permutation_conjugate_",ToString[n],"(
const double *__restrict__ x, const double *__restrict__ p,
double *__restrict__ y) {\n",MapThread[StringJoin,{
"y["<>ToString[#]<>"] = "&/@Range[0,n (n+1)/2-1],
toCCode@*Expand/@Join@@MapThread[Take,
{p.packedMatrix["x",n].Transpose[p],Range[n]}],
ConstantArray[";\n",n (n+1)/2]}],"}\n"]];

StringJoin@Riffle[Table[StringJoin[
determinantInverseFunctionCCode[n],"\n\n",
kineticTraceFunctionCCode[n],"\n\n",
quadraticFormFunctionCCode[n],"\n\n",
permutationConjugateFunctionCCode[n]],
{n,7}],"\n\n"]//CopyToClipboard;

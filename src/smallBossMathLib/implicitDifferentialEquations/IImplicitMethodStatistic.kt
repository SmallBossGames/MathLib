package smallBossMathLib.implicitDifferentialEquations

interface IImplicitMethodStatistic {
    val stepsCount: Int
    val evaluationsCount: Int
    val jacobiEvaluationsCount: Int
    val returnsCount: Int
}
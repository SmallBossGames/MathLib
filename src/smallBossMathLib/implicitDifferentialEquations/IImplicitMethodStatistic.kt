package smallBossMathLib.implicitDifferentialEquations

interface IImplicitMethodStatistic {
    val isLowStepSizeReached:Boolean
    val isHighStepSizeReached: Boolean
    val stepsCount: Int
    val evaluationsCount: Int
    val jacobiEvaluationsCount: Int
    val returnsCount: Int
}
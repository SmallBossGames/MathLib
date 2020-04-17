package smallBossMathLib.explicitDifferentialEquations

interface IExplicitMethodStepInfo {
    val isLowStepSizeReached:Boolean
    val isHighStepSizeReached: Boolean
    val stepsCount: Int
    val evaluationsCount: Int
    val returnsCount: Int
}
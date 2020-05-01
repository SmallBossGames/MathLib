package smallBossMathLib.explicitDifferentialEquations

interface IExplicitMethodStepState {
    val isLowStepSizeReached:Boolean
    val isHighStepSizeReached: Boolean
}
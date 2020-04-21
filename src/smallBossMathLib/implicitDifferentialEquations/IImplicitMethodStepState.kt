package smallBossMathLib.implicitDifferentialEquations

interface IImplicitMethodStepState {
    val maxEigenvalue: Double
    val minEigenvalue: Double
    val isLowStepSizeReached: Boolean
    val isHighStepSizeReached: Boolean
    val freezeJacobiStepsCount: Int
}
package smallBossMathLib.implicitDifferentialEquations

class ImplicitMethodStepState(
    override val isLowStepSizeReached: Boolean,
    override val isHighStepSizeReached: Boolean,
    override val freezeJacobiStepsCount: Int
) : IImplicitMethodStepState {
    override val maxEigenvalue: Double
        get() = TODO("Not yet implemented")
    override val minEigenvalue: Double
        get() = TODO("Not yet implemented")

}
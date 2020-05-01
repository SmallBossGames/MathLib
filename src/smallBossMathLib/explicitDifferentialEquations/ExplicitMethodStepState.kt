package smallBossMathLib.explicitDifferentialEquations

data class ExplicitMethodStepState(
    override var isLowStepSizeReached: Boolean,
    override var isHighStepSizeReached: Boolean
) : IExplicitMethodStepState
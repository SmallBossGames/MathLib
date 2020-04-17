package smallBossMathLib.explicitDifferentialEquations

import smallBossMathLib.explicitDifferentialEquations.IExplicitMethodStepInfo

data class ExplicitMethodStepInfo(
    override var isLowStepSizeReached: Boolean,
    override var isHighStepSizeReached: Boolean,
    override var stepsCount: Int,
    override var evaluationsCount: Int,
    override var returnsCount: Int
) : IExplicitMethodStepInfo
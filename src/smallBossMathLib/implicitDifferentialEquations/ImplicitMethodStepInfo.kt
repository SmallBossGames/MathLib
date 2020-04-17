package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.implicitDifferentialEquations.IImplicitMethodStepInfo

data class ImplicitMethodStepInfo(
    override var isLowStepSizeReached: Boolean,
    override var isHighStepSizeReached: Boolean,
    override var stepsCount: Int,
    override var evaluationsCount: Int,
    override var jacobiEvaluationsCount: Int,
    override var returnsCount: Int
) : IImplicitMethodStepInfo
package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.shared.Matrix2D

data class ImplicitMethodStepState(
    override val jacobiMatrix: Matrix2D,
    override var isLowStepSizeReached: Boolean,
    override var isHighStepSizeReached: Boolean,
    override var freezeJacobiStepsCount: Int
) : IImplicitMethodStepState
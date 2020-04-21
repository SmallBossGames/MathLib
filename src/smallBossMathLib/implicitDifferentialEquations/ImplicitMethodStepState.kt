package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.shared.Matrix2D

data class ImplicitMethodStepState(
    override var isLowStepSizeReached: Boolean,
    override var isHighStepSizeReached: Boolean,
    override var maxEigenvalue: Double,
    override var minEigenvalue: Double,
    override var freezeJacobiStepsCount: Int
) : IImplicitMethodStepState
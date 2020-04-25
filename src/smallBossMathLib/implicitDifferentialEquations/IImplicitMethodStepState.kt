package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.shared.Matrix2D

interface IImplicitMethodStepState {
    val jacobiMatrix: Matrix2D
    val isLowStepSizeReached: Boolean
    val isHighStepSizeReached: Boolean
    val freezeJacobiStepsCount: Int
}
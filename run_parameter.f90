program main
    use BuFlowModule  
!    use meshdeformationn  
    implicit none
    ! 定义变形参数（1×4矩阵）
!   real(kind=8) :: data_4D137(1,4)
    
    ! 声明数据结构时，是(1,4)的矩阵，是二维数组。下面的reshape就是把一维数组[-0.75d0, 0.485d0, -0.588d0, -0.912d0]重塑为二维数组（第一维度是1）
!    data_4D137 = reshape([0.0d0, 0.0d0, 0.0d0, 0.0d0], [1,4])
    ! 可选参数值（注释保留原格式）
    ! data_4D137 = reshape([0.0d0, 0.0d0, 0.0d0, 0.0d0], [1,4])
    ! data_4D137 = reshape([-0.221d0, -0.382d0, 0.941d0, -0.456d0], [1,4])

    ! 网格变形
!    call airfoil_deformation_HH(data_4D137)
    !计算CFD
    call compute_CFD_RANS
end program main




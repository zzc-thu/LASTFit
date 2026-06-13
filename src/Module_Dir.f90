module directory_utils
    use, intrinsic :: iso_c_binding
    implicit none

    ! Bind the mkdir function to the C
    interface
        function c_mkdir(path, mode) bind(C, name="mkdir")
            use iso_c_binding
            implicit none
            integer(c_int) :: c_mkdir
            character(kind=c_char), intent(in) :: path(*) ! C风格字符串（Null结尾）
            integer(c_int), intent(in) :: mode           ! 权限模式
        end function c_mkdir
    end interface

contains

    ! 子程序：检查并创建目录
    subroutine ensure_directory_exists(dir_name, status)
        ! 输入参数
        character(len=*), intent(in) :: dir_name  ! 输入目录名称
        ! 输出参数
        integer, intent(out) :: status           ! 返回状态：0 成功，非 0 失败

        logical :: dir_exists                    ! 用于 `inquire` 检查目录是否存在
        integer(c_int) :: mkdir_status           ! mkdir 执行的返回状态
        integer(c_int), parameter :: mode = int(o'755', c_int)  ! 默认权限

        ! 确保目录名以 C 字符串格式传入（带 Null Terminator）
        character(len=len(dir_name)+1) :: c_dir_name
        c_dir_name = trim(dir_name) // c_null_char

        ! 检查目录是否已存在
        inquire(file=dir_name, exist=dir_exists)
        if (dir_exists) then
            status = 0  ! 如果目录已经存在，直接返回成功
            return
        end if

        ! 如果目录不存在，则创建
        mkdir_status = c_mkdir(c_dir_name, mode)
        if (mkdir_status == 0) then
            status = 0  ! 目录创建成功
        else
            status = 1  ! 目录创建失败
        end if
    end subroutine ensure_directory_exists

end module directory_utils

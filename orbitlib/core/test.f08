program sample0
    use orbit_func
    use orbit_base
    use orbit_debri
    implicit none

    type(TOrbit) :: orbit

    ! モジュール変数DEBRISにデータを読み込む
    call load_debri(num_debris=1, file_tle="data/debri_elements.txt", file_rcs="data/RCS_list.txt", check_size=.false.)

    orbit = DEBRIS(1)%orbit

    print *, "size=", size(DEBRIS)
    print *, "a=", orbit%semimajor_axis
end program sample0

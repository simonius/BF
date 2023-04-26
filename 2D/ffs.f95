! This file is the main driver routine for the forward facing step problem
! The results are saved in a single HDF5 File


program ffs
        use consts
        use fluxes
        use solvers
        use HDF5
        use FFSdomain

        implicit none

!       IO Parameters and ids
!       Filename under which the result will be saved
        Character(Len=30) :: filename = "ffsout"
!       String that descrives the grid used
        Character(Len=30) :: gridstring
!       Name of the dataset into which the array is saved
        Character(Len=20) :: dsetname = "u"
!       File Identifier, Group identidier, Dataset identifier, Dataspace Identifier
        Integer(HSIZE_T) :: file_id, group_id, dset_id, dspace_id
!       takes returned errors
        Integer :: error

!       Dataset information
        Integer(HSIZE_T), Dimension(4) :: dims
        Integer :: rank = 4

! solver matrices
        double precision, dimension(nx, ny, ncons) :: u0ar
        double precision, dimension(nsnaps, nx, ny, ncons) :: uar
        double precision :: tmax = 4.0
        double precision :: tsnapdist = 0.0

!
        integer :: i, j, k

!       writing grid size into the grid string
        write(gridstring, *) nx

!       combining file name, grid string and file ending
        filename = trim(filename) // trim(adjustl(gridstring)) // ".h5"

!       writing out the filename for verification
        write (*,*) filename

!       Initialize the dimensions of the HDF5 file
        dims(1) = nsnaps
        dims(2) = nx
        dims(3) = ny
        dims(4) = ncons

        tsnapdist = tmax/dble(nsnaps)

!       Open the HDF5 Fortran interfce
        call h5open_f(error)

!       Create an HDF5 File
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!       Create the Dataspace in the HDF5 File
        call h5screate_simple_f(rank,  dims, dspace_id, error)

!       Create the dataset inside this dataspace
        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

!       Initializing
        i = init()

!       set boundary implementation pointers
        sanitizeux=>sanitizeffsx
        sanitizeuy=>sanitizeffsy

!       set distance implementations
        updist=>updistfunc
        downdist=>downdistfunc

!       set initial conditions.
        call initffs(u0ar)
        
!       solving
        call solEGT(uar(1, :, :, :), u0ar, tsnapdist)
        do i=2, nsnaps
                call solEGT(uar(i, :, :, :), uar(i-1, :, :, :), tsnapdist)
        end do

!       saving the result in the HDF5 file
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, uar, dims, error)

        ! closing all HDF5 things
        CALL h5dclose_f(dset_id, error)
        CALL h5sclose_f(dspace_id, error)
        CALL h5fclose_f(file_id, error)
        CALL h5close_f(error)


end program ffs


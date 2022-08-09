subroutine getBasisPtSurface(u, v, tu, tv, ku, kv, vals, col_ind, istart, l_index, nctlu, nctlv, nnz)
    use precision
    implicit none

    ! Input/Output
    integer, intent(in) :: ku, kv, nctlu, nctlv, nnz, istart
    real(kind=realType), intent(in) :: u, v
    real(kind=realType), intent(in) :: tu(nctlu + ku), tv(nctlv + kv)
    integer, intent(in) :: l_index(Nctlv, Nctlu)
    real(kind=realType), intent(inout) :: vals(nnz)
    integer, intent(inout) :: col_ind(nnz)

    ! Working
    real(kind=realType) :: basisu(ku), basisv(kv)
    integer :: ileftu, ileftv
    integer :: ii, jj, counter, start

    ! Get u interval
    call findSpan(u, ku, tu, nctlu, ileftu)
    call basis(tu, nctlu, ku, u, ileftu, basisu)

    ! Get v interval
    call findSpan(v, kv, tv, nctlv, ileftv)
    call basis(tv, nctlv, kv, v, ileftv, basisv)

    counter = 0
    do ii = 1, ku
        do jj = 1, kv
            ! Get the local row/col for this surface
            start = istart + counter + 1
            col_ind(start) = l_index(ileftv - kv + jj, ileftu - ku + ii)
            vals(start) = basisu(ii) * basisv(jj)
            counter = counter + 1
        end do
    end do
end subroutine getBasisPtSurface

subroutine getBasisPtVolume(u, v, w, tu, tv, tw, ku, kv, kw, vals, col_ind, istart, l_index, nctlu, nctlv, nctlw, nnz)
    use precision
    implicit none

    ! Input
    integer, intent(in) :: ku, kv, kw, nctlu, nctlv, nctlw, nnz, istart
    real(kind=realType), intent(in) :: u, v, w
    real(kind=realType), intent(in) :: tu(nctlu + ku), tv(nctlv + kv), tw(nctlw + kw)
    integer, intent(in) :: l_index(Nctlw, Nctlv, Nctlu)

    ! Output
    real(kind=realType), intent(inout) :: vals(nnz)
    integer, intent(inout) :: col_ind(nnz)

    ! Working
    real(kind=realType) :: basisu(ku), basisv(kv), basisw(kw)
    integer :: ileftu, ileftv, ileftw
    integer :: ii, jj, kk, counter, start

    ! Get u interval
    call findSpan(u, ku, tu, nctlu, ileftu)
    call basis(tu, nctlu, ku, u, ileftu, basisu)

    ! Get v interval
    call findSpan(v, kv, tv, nctlv, ileftv)
    call basis(tv, nctlv, kv, v, ileftv, basisv)

    ! Get w interval
    call findSpan(w, kw, tw, nctlw, ileftw)
    call basis(tw, nctlw, kw, w, ileftw, basisw)

    counter = 0
    do ii = 1, ku
        do jj = 1, kv
            do kk = 1, kw
                ! Get the local row/col for this surface
                start = istart + counter + 1
                col_ind(start) = l_index(ileftw - kw + kk, ileftv - kv + jj, ileftu - ku + ii)

                vals(start) = basisu(ii) * basisv(jj) * basisw(kk)
                counter = counter + 1
            end do
        end do
    end do
end subroutine getBasisPtVolume


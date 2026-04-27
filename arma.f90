module arma_mod
use kind_mod, only: dp
implicit none
private
public :: ar_fit,ma1_acf,ma2_acf,ma_acf,ma_rmse,ma_pred,ar_pred,ar_coeff,ar_err, &
          arma_pred,arma_rmse,arma_rmse_i1_i2,arma_psi,arma_acv,arma_acf,arma_1_1_acf, &
          ma_2_acf,ma_acv,arma_coeff_string,ar_2_acf,print_arma,print_ar_fit
character (len=*), parameter :: mod_str = "in arma_mod::"
logical, parameter :: bound_arma_pred = .true.
real(kind=dp), parameter :: arma_pred_min = -10.0_dp, arma_pred_max = 10.0_dp, bad_real = -999.0_dp
integer, parameter :: istdout = 6
interface default
   module procedure default_integer,default_logical
end interface default
interface set_alloc
   module procedure set_alloc_real_vec
end interface set_alloc
contains
function ar_err(xx,nar) result(xerr)
! return the errors from fitting an autoregressive model of order nar to xx(:)
real(kind=dp), intent(in) :: xx(:)
integer      , intent(in) :: nar
real(kind=dp)             :: xerr(size(xx))
real(kind=dp)             :: xar(nar)
xar = ar_coeff(xx,nar)
xerr = xx - ar_pred(xx,xar)
end function ar_err
!
function ar_pred(xx,xar) result(yy)
! return predictions from a zero-mean autoregressive (AR) model with AR coefficients xar(:)
real(kind=dp), intent(in) :: xx(:)
real(kind=dp), intent(in) :: xar(:)
real(kind=dp)             :: yy(size(xx))
integer                   :: i,nw,n
n  = size(xx)
nw = size(xar)
if (n < 1) return
if (size(yy) /= n) then
   yy = huge(xx)
   return
else if (nw < 1) then
   yy = 0.0_dp
   return
end if
yy(1) = 0.0_dp
do i=2,n
   if (i > nw) then
      yy(i) = sum(xx(i-1:i-nw:-1)*xar)
   else
      yy(i) = sum(xx(i-1:1:-1)*xar(:i-1))
   end if
end do
end function ar_pred
!
function ar_coeff(xx,nar) result(xar)
! return autoregressive coefficients calculated from call to ar_fit
real(kind=dp), intent(in) :: xx(:)
integer      , intent(in) :: nar
real(kind=dp)             :: xar(nar)
real(kind=dp)             :: xms
call ar_fit(xx,xms,xar)
end function ar_coeff
!
subroutine print_ar_fit(xx,nar_max,outu,title,fmt_header,best_order,best_ar, &
                        print_each_fit)
! fit AR (autoregressive) models to xx(:), up to order nar_max
real(kind=dp)    , intent(in)            :: xx(:)
integer          , intent(in) , optional :: nar_max
integer          , intent(in) , optional :: outu
character (len=*), intent(in) , optional :: title,fmt_header
integer          , intent(out), optional :: best_order
real(kind=dp)    , intent(out), optional, allocatable :: best_ar(:)
logical          , intent(in) , optional :: print_each_fit
integer                                  :: nar_max_,iar,n,outu_,best_order_
real(kind=dp)                            :: xms,rmse,aicc,aicc_lowest
real(kind=dp)    , allocatable           :: xar(:)
logical                                  :: print_each_fit_,print_best_fit_
print_each_fit_ = default(.true.,print_each_fit)
print_best_fit_ = .false.
nar_max_ = default(min(size(xx)/4,10),nar_max)
outu_ = default(istdout,outu)
n = size(xx)
allocate (xar(nar_max_))
if (present(fmt_header)) write (outu_,fmt_header)
if (present(title)) write (outu_,"(a)") trim(title)
best_order_ = 0
if (present(best_ar)) allocate (best_ar(0))
if (print_each_fit_) &
   write (outu_,"(/,3a8,100(:,2x,'AR(',i2.2,')'))") "order","AICC","RMSE",(iar,iar=1,nar_max_)
do iar=0,nar_max_
   call ar_fit(xx,xms,xar(:iar))
   rmse = sqrt(xms)
   aicc = 2*log(rmse)+(n+iar)/(n-iar-2.0_dp)
   if (iar == 0) then
      aicc_lowest = aicc
   else if (aicc < aicc_lowest) then
      best_order_ = iar
      aicc_lowest = aicc
      if (present(best_ar)) call set_alloc(xar(:iar),best_ar)
   end if
   if (print_each_fit_) write (outu_,"(i8,100f8.4)") iar,aicc,rmse,xar(:iar)
end do
if (print_best_fit_) write (outu_,"('best order = ',i0)") best_order_
if (present(best_order)) best_order = best_order_
end subroutine print_ar_fit
!
subroutine ar_fit(data,xms,xar,xpacf)
! fit an autoregressive model using the Burg algorithm
real(dp), intent(in)  :: data(:)
real(dp), intent(out) :: xms       ! mean_squared error
real(dp), intent(out) :: xar(:)    ! AR coefficients
real(dp), intent(out), optional :: xpacf(:) ! partial autocorrelations
integer                         :: k,m,n
real(dp)                        :: denom,pneum
real(dp), dimension(size(data)) :: wk1,wk2,wktmp
real(dp), dimension(size(xar))  :: wkm
logical , parameter             :: print_pacf = .false.
character (len=*), parameter    :: msg = mod_str // "ar_fit, "
if (present(xpacf)) then
   if (size(xpacf) /= size(xar)) then
      write (*,*) msg,"size(xpacf), size(xar) =",size(xpacf),size(xar)," should be equal, STOPPING"
      stop
   end if
end if
m          = size(xar)
n          = size(data)
xms        = dot_product(data,data)/n
! print*,"in memcof, m, n=",m,n ! debug
if (m < 1) return
wk1(1:n-1) = data(1:n-1)
wk2(1:n-1) = data(2:n)
do k = 1,m
   pneum      = dot_product(wk1(1:n-k),wk2(1:n-k))
   denom      = dot_product(wk1(1:n-k),wk1(1:n-k)) + dot_product(wk2(1:n-k),wk2(1:n-k))
   xar(k)     = 2.0_dp*pneum/denom
   xms        = xms*(1.0_dp-xar(k)**2)
   xar(1:k-1) = wkm(1:k-1)-xar(k)*wkm(k-1:1:-1)
   if (present(xpacf)) xpacf(k) = xar(k)
   if (print_pacf) then
      write (*,"(1x,'k, ar(k) =',1x,i6,1x,100f9.4)") k,xar(1:k)
   end if
   if (k == m) return
   wkm(1:k)     = xar(1:k)
   wktmp(2:n-k) = wk1(2:n-k)
   wk1(1:n-k-1) = wk1(1:n-k-1) - wkm(k)*wk2(1:n-k-1)
   wk2(1:n-k-1) = wk2(2:n-k)   - wkm(k)*wktmp(2:n-k)
end do
stop "should not get here in memcof"
end subroutine ar_fit
!
elemental function ma1_acf(ma1) result(acf1)
! return ACF(1) of an MA(1) process using formula from https://onlinecourses.science.psu.edu/stat510/node/48/
real(kind=dp), intent(in) :: ma1
real(kind=dp)             :: acf1
acf1 = ma1/(1+ma1**2)
end function ma1_acf
!
pure function ma2_acf(ma1,ma2) result(xacf)
! return ACF(1:2) of an MA(2) process using formula from https://onlinecourses.science.psu.edu/stat510/node/48/
real(kind=dp), intent(in) :: ma1,ma2
real(kind=dp)             :: xacf(2)
xacf = [ma1 + ma1*ma2, ma2]/(1 + ma1**2 + ma2**2)
end function ma2_acf
!
pure function ma_acf(ma) result(xacf)
! return ACF of an MA process using formula from https://www.stat.ncsu.edu/people/bloomfield/courses/ST730/slides/SnS-03-3.pdf 
real(kind=dp), intent(in) :: ma(:)
real(kind=dp)             :: xacf(size(ma))
real(kind=dp)             :: denom
integer                   :: ilag,j,nma
! j = 1
nma   = size(ma)
denom = 1.0_dp + sum(ma**2)
do ilag=1,nma
   xacf(ilag) = (ma(ilag) + sum([(ma(j)*ma(ilag+j),j=1,nma-ilag)]))/denom
end do
end function ma_acf
!
function ma_rmse(xx,ma) result(rms_err)
! compute root-mean-squared error of moving average model with coefficients ma(:) applied to time series xx(:)
real(kind=dp), intent(in) :: xx(:),ma(:)
real(kind=dp)             :: rms_err
rms_err = rms(xx-ma_pred(xx,ma))
end function ma_rmse
!
function arma_rmse(xx,ar,ma) result(rms_err)
! compute root-mean-squared error of ARMA model with autoregressive coefficients ar(:) and moving average coefficients ma(:) applied to time series xx(:)
real(kind=dp), intent(in) :: xx(:),ar(:),ma(:)
real(kind=dp)             :: rms_err
! print*,"minval(arma_pred), maxval(arma_pred)=",minval(arma_pred(xx,ar,ma)),maxval(arma_pred(xx,ar,ma))
rms_err = rms(xx-arma_pred(xx,ar,ma))
end function arma_rmse
!
function arma_rmse_i1_i2(xx,ar,ma,i1,i2) result(rms_err)
! compute root-mean-squared error of ARMA model with autoregressive coefficients ar(:) and moving average coefficients ma(:) applied to time series xx(:)
real(kind=dp), intent(in) :: xx(:),ar(:),ma(:)
integer      , intent(in), optional :: i1,i2
real(kind=dp)             :: rms_err
real(kind=dp)             :: xerr(size(xx))
integer                   :: n,i1_,i2_
n       = size(xx)
xerr    = xx-arma_pred(xx,ar,ma)
if (present(i1)) i1_ = max(1,i1)
if (present(i2)) i2_ = min(n,i2)
if (present(i1) .and. present(i2)) then
   rms_err = rms(xerr(i1_:i2_))
else if (present(i1)) then
   rms_err = rms(xerr(i1_:))
else if (present(i2)) then
   rms_err = rms(xerr(:i2_))
else
   rms_err = rms(xerr)
end if
end function arma_rmse_i1_i2
!
function ma_pred(xx,ma) result(xpred)
! return predictions of moving average model
real(kind=dp), intent(in) :: xx(:),ma(:)
real(kind=dp)             :: xpred(size(xx))
integer                   :: i,n,nma
real(kind=dp)             :: ee(1-size(ma):size(xx))
n   = size(xx)
if (n < 1) return
nma = size(ma)
ee  = 0.0_dp
do i=1,n
   xpred(i) = sum(ma(nma:1:-1)*ee(i-nma:i-1))
   ee(i)    = xx(i) - xpred(i)
end do
end function ma_pred
!
function arma_pred(xx,ar,ma) result(xpred)
! return predictions of autoregressive moving average model
real(kind=dp), intent(in) :: xx(:),ar(:),ma(:)
real(kind=dp)             :: xpred(size(xx))
integer                   :: i,n,nar,nar_use,nma
real(kind=dp)             :: ee(1-size(ma):size(xx))
n   = size(xx)
if (n < 1) return
nar = size(ar)
nma = size(ma)
ee  = 0.0_dp
do i=1,n
   nar_use  = min(nar,i-1)
   xpred(i) = sum(ma(nma:1:-1)*ee(i-nma:i-1)) + sum(ar(:nar_use)*xx(i-1:i-nar_use:-1))
   if (bound_arma_pred) xpred(i) = min(arma_pred_max,max(arma_pred_min,xpred(i)))
   ee(i)    = xx(i) - xpred(i)
end do
end function arma_pred
!
pure function rms(xx) result(xrms)
! compute the root-mean-squared value of xx
real(kind=dp), intent(in) :: xx(:)
real(kind=dp)             :: xrms
integer                   :: n
n = size(xx)
if (n > 0) then
   xrms = sqrt(sum(xx**2)/n)
else
   xrms = -1.0
end if 
end function rms
!
function arma_psi(ar,ma,nlags) result(psi)
! compute the psi (moving average) weights corresponding to an ARMA process
real(kind=dp), intent(in) :: ar(:),ma(:) ! autoregressive and moving average weights
integer      , intent(in) :: nlags       ! # of psi weights to compute
real(kind=dp)             :: psi(nlags)
integer                   :: j,k,p,q
if (nlags < 1) return
psi = 0.0_dp
p   = size(ar)
q   = size(ma)
do j=1,nlags ! min(nlags,max(p,q+1))
!   print*,"j=",j
   psi(j) = 0.0_dp
   if (j <= q) psi(j) = psi(j) + ma(j)
   if (j <= p) psi(j) = psi(j) + ar(j) ! + sum(ar(1:j-1)*psi(
   do k=1,min(j-1,p)
      psi(j) = psi(j) + ar(k)*psi(j-k)
   end do
end do
end function arma_psi
!
function arma_acv(ar,ma,nlags) result(acv)
! compute the autocovariance corresponding to an ARMA process
real(kind=dp), intent(in) :: ar(:),ma(:) ! autoregressive and moving average weights
integer      , intent(in) :: nlags       ! # of psi weights to compute
real(kind=dp)             :: acv(0:nlags)
real(kind=dp)             :: psi(0:size(ma))
integer                   :: j,h,p,q
p   = size(ar)
q   = size(ma)
psi = [0.0_dp,arma_psi(ar,ma,q)]
do h=0,nlags
   if (h == 0) then
      acv(h) = 1.0_dp
   else if (h <= q) then
      acv(h) = acv(0)*psi(h)
   else
      acv(h) = 0.0_dp
   end if
   if (h < max(p,q+1)) then
      do j=1,min(p,h)
         acv(h) = acv(h) + ar(j)*acv(h-j)
      end do
      do j=max(1,h),q
         acv(h) = acv(h) + ma(j)*psi(j-h)
      end do
   else
      do j=1,p
         acv(h) = acv(h) + ar(j)*acv(h-j)
      end do
   end if
end do
end function arma_acv
!
function arma_acf(ar,ma,nlags) result(acf)
! compute the autocovariance corresponding to an ARMA process
real(kind=dp), intent(in) :: ar(:),ma(:) ! autoregressive and moving average weights
integer      , intent(in) :: nlags       ! # of ACF to compute
real(kind=dp)             :: acf(nlags)
real(kind=dp)             :: acv(0:nlags)
integer      , parameter  :: nlags_psi = 50
real(kind=dp)             :: psi(nlags_psi),psi_acf(nlags_psi)
if (.false.) then
   acv = arma_acv(ar,ma,nlags)
   acf = acv(1:nlags)/acv(0)
end if
psi = arma_psi(ar,ma,nlags_psi)
psi_acf = ma_acf(psi)
if (nlags < nlags_psi) then
   acf = psi_acf(:nlags)
else
   acf(:nlags_psi) = psi_acf
   acf(nlags_psi+1:) = 0.0_dp
end if  
end function arma_acf
!
function arma_1_1_acf(ar,ma,nlags) result(acf)
! return ACF of an ARMA(1,1) process
! Shumway and Stoffer (2017), p96, example 3.14
real(kind=dp), intent(in) :: ar,ma  ! autoregressive and moving average weights
integer      , intent(in) :: nlags  ! # of ACF to compute
real(kind=dp)             :: acf(nlags)
real(kind=dp)             :: denom
integer                   :: i
if (nlags < 1) return
denom = 1.0_dp + 2*ar*ma + ma**2
if (denom /= 0.0_dp) then
   acf(1) = ((ar + ma)*(1.0_dp + ar*ma))/denom
else
   acf(1) = bad_real
end if
do i=2,nlags
   acf(i) = ar*acf(i-1)
end do
end function arma_1_1_acf
!
function ar_2_acf(ar1,ar2,nlags) result(acf)
! return the ACF of an AR(2) process
real(kind=dp), intent(in) :: ar1,ar2  ! moving average coefficients
integer      , intent(in) :: nlags  ! # of ACF to compute
real(kind=dp)             :: acf(nlags)
real(kind=dp)             :: denom
integer                   :: i
if (nlags < 1) return
acf    = 0.0_dp
denom  = 1.0_dp - ar2
if (denom == 0.0_dp) then
   acf = bad_real
   return
end if
acf(1) = ar1/denom
if (nlags < 2) return
acf(2) = ar1**2/denom + ar2
do i=3,nlags
   acf(i) = ar1*acf(i-1) + ar2*acf(i-2)
end do
end function ar_2_acf
!
function ma_2_acf(ma1,ma2,nlags) result(acf)
! return the ACF of an MA(2) process
! Applied Time Series Analysis, 2.1 Moving Average Models (MA models), https://onlinecourses.science.psu.edu/stat510/node/48/
real(kind=dp), intent(in) :: ma1,ma2  ! moving average coefficients
integer      , intent(in) :: nlags  ! # of ACF to compute
real(kind=dp)             :: acf(nlags)
real(kind=dp)             :: denom
if (nlags < 1) return
acf    = 0.0_dp
denom  = 1.0_dp + ma1**2 + ma2**2
acf(1) = (ma1 + ma1*ma2)/denom
if (nlags > 1) acf(2) = ma2/denom
end function ma_2_acf
!
function ma_acv(ma,nlags) result(acv)
! Shumway and Stoffer (2017) (3.42) p94
! compute the autocovariance corresponding to an MA process, including the 0th lag ACV
real(kind=dp), intent(in) :: ma(:) ! moving average weights
integer      , intent(in) :: nlags ! # of psi weights to compute
real(kind=dp)             :: acv(0:nlags)
real(kind=dp)             :: xma(0:size(ma))
integer                   :: j,h,q
xma = [1.0_dp,ma]
q   = size(ma)
acv = 0.0_dp
do h=0,min(q,nlags)
   do j=0,q-h
      acv(h) = acv(h) + xma(j)*xma(j+h)
   end do
end do
end function ma_acv
!
pure function arma_coeff_string(ar,ma) result(text)
real(kind=dp), intent(in) :: ar(:),ma(:)
character (len=1000) :: text
write (text,"('AR=[',f0.3,100(1x,f0.3))") ar
write (text,"(a,'] MA=[',f0.3,100(1x,f0.3))") trim(text),ma
text = trim(text) // "]"
end function arma_coeff_string
!
subroutine print_arma(ar,ma,outu)
real(kind=dp), intent(in) :: ar(:),ma(:)
integer      , intent(in), optional :: outu
integer                             :: outu_
if (present(outu)) then
   outu_ = outu
else
   outu_ = istdout
end if
if (size(ar) > 0) write (*,"(  'AR =',100f8.4)") ar
if (size(ma) > 0) write (*,"(  'MA =',100f8.4)") ma
end subroutine print_arma
!
elemental function default_integer(def,opt) result(ii)
! return opt if it is present, otherwise def (integer arguments and result)
integer, intent(in)           :: def
integer, intent(in), optional :: opt
integer                       :: ii
if (present(opt)) then
   ii = opt
else
   ii = def
end if
end function default_integer
!
elemental function default_logical(def,opt) result(tf)
! return opt if it is present, otherwise def (logical arguments and result)
logical, intent(in)           :: def
logical, intent(in), optional :: opt
logical                       :: tf
if (present(opt)) then
   tf = opt
else
   tf = def
end if
end function default_logical
!
subroutine set_alloc_real_vec(xx,yy,nsize)
! for real vectors xx(:) and yy(:), allocate yy and set it to xx
real(kind=dp), intent(in)               :: xx(:)
real(kind=dp), intent(out), allocatable :: yy(:)
integer      , intent(out), optional    :: nsize
integer                                 :: ni
ni = size(xx)
allocate (yy(ni))
yy = xx
if (present(nsize)) nsize = ni
end subroutine set_alloc_real_vec
!
end module arma_mod

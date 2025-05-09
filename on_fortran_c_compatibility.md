# Gotcha #1

You need to have: `real(c_double), intent(in) :: selfenergy(*)`.
- `real(c_double)`, with `c_double` is mandatory
- `intent(in)` is *not* mandatory
- `selfenergy(*)` notice the `(*)` is mandatory. A `selfenergy(:)` or `selfenergy(:,:)` will not work.
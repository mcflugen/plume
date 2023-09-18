import numpy as np
cimport numpy as np
cimport cython

from cython_gsl cimport *
from libc.stdlib cimport malloc, free

# cdef extern from "stdlib.h":
#     void *malloc(size_t size)
#     int free(void*)
#     int sizeof()

ctypedef double * double_ptr
ctypedef void * void_ptr


cdef double f(double x, void * params) noexcept nogil:
    cdef double x0 = (<double_ptr> params)[0]
    cdef double y0 = (<double_ptr> params)[1]
    cdef double w = (<double_ptr> params)[2]
    cdef double c = (<double_ptr> params)[3]
    cdef double n = (<double_ptr> params)[4]
    cdef double denom = (pow(w * pow(x / (c * w), 1. / n) - y0, 2.) +
                         pow(x - x0, 2.))

    return ((pow(x / (c * w), 1. / n - 1.) *
             (w * pow(x / (c * w), 1. / n) - y0) /
             (n * c) + (x - x0)) / sqrt(denom))


cdef double r(double x, double x0, double y0, double w, double c,
              double n) nogil:
    return sqrt(pow(x - x0, 2.) + pow(centerline_y(x, w, c, n) - y0,  2.))


cdef void find_bounds(void * params, double *bounds) nogil:
    cdef double x0 = (<double_ptr> params)[0]
    cdef double y0 = (<double_ptr> params)[1]
    cdef double w = (<double_ptr> params)[2]
    cdef double c = (<double_ptr> params)[3]
    cdef double n = (<double_ptr> params)[4]
    cdef double lower, upper, tmp
    cdef int n_samples = 10
    cdef double x_, dx, x_of_min
    cdef double distance, min_distance

    lower = x0
    if lower < 0.:
        lower = 0.

    if y0 < 0:
        upper = 0.
    else:
        upper = centerline_x(y0, w, c, n)

    if lower > upper:
        tmp = lower
        lower = upper
        upper = tmp

    upper *= 1.2
    lower *= .8

    dx = (upper - lower) / (n_samples - 1)

    x_ = lower
    min_distance = r(x_, x0, y0, w, c, n)
    x_of_min = x_
    for _ in range(1, n_samples):
        x_ += dx
        distance = r(x_, x0, y0, w, c, n)
        if distance < min_distance:
            min_distance = distance
            x_of_min = x_

    bounds[0] = x_of_min - dx
    bounds[1] = x_of_min + dx

    if bounds[0] < 0.:
        bounds[0] = 0.


@cython.cdivision(True)
cdef double centerline_y(double x, double rw, double c, double n) nogil:
    return rw * pow(x / (rw * c), 1. / n)


@cython.cdivision(True)
cdef double centerline_x(double y, double rw, double c, double n) nogil:
    return rw * c * pow(y / rw, n)


@cython.cdivision(True)
cdef double g(double x, void * params) noexcept nogil:
    cdef double w = (<double_ptr> params)[0]
    cdef double c = (<double_ptr> params)[1]
    cdef double n = (<double_ptr> params)[2]
    cdef double m = (<double_ptr> params)[3]
    cdef double b = (<double_ptr> params)[4]

    return (w * pow(x / (w * c), 1. / n)) - (m * x + b)


cdef void c_y_rotated(double * xc, const int n_points, const double rw,
                      const double c, const double n, const double a,
                      double * yc) nogil:
    cdef int GSL_CONTINUE = -2
    cdef int MAX_ITERS = 100
    cdef double SIN_A = sin(a)
    cdef double COS_A = cos(a)
    cdef gsl_root_fsolver *s
    cdef gsl_root_fsolver_type *T
    cdef gsl_function F
    cdef double params[5]
    cdef double bounds[2]
    cdef int status
    cdef int iters
    cdef int point
    cdef double x1, y1
    cdef double m = tan(M_PI * .5 - a)
    cdef double b

    if fabs(SIN_A) < 1e-12:
        for point in range(n_points):
            yc[point] = centerline_y(COS_A * xc[point], rw, c, n) * COS_A
        return

    if fabs(COS_A) < 1e-12:
        for point in range(n_points):
            yc[point] = centerline_x(- SIN_A * xc[point], rw, c, n) * SIN_A
        return

    params[0] = rw
    params[1] = c
    params[2] = n
    params[3] = m

    bounds[0] = 0.

    F.function = &g
    F.params = params

    T = gsl_root_fsolver_brent
    s = gsl_root_fsolver_alloc(T)

    for point in range(n_points):
        x1 = xc[point] * COS_A
        y1 = - xc[point] * SIN_A
        b = y1 - m * x1

        params[4] = b

        bounds[1] = - b / m

        gsl_root_fsolver_set(s, &F, bounds[0], bounds[1])

        status = GSL_CONTINUE
        iters = 0
        while status == GSL_CONTINUE:
            status = gsl_root_fsolver_iterate(s)
            x_intercept = gsl_root_fsolver_root(s)

            bounds[0] = gsl_root_fsolver_x_lower(s)
            bounds[1] = gsl_root_fsolver_x_upper(s)

            status = gsl_root_test_interval(bounds[0], bounds[1], 0., 0.001)

            iters += 1
            if iters > MAX_ITERS:
                break

        y_intercept = centerline_y(x_intercept, rw, c, n)
        yc[point] = (y_intercept - y1) / COS_A

    gsl_root_fsolver_free(s)

    return


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def y_rotated(np.ndarray[np.float_t, ndim=1, mode="c"] xc not None,
              double rw, double c, double n, double a,
              np.ndarray[np.float_t, ndim=1, mode="c"] yc not None):
    c_y_rotated(&xc[0], len(xc), rw, c, n, a, &yc[0])


cdef int nearest_is_lower(void * params, double * bounds) nogil:
    f_lower = f(bounds[0], params)
    f_upper = f(bounds[1], params)

    return f_lower * f_upper > 0


cdef void c_path_lengths(double * bounds, const int n_points, const double rw,
                         const double c, const double n, const double a,
                         double * lengths) nogil:
    cdef int N_SAMPLES = 100
    cdef int point
    cdef double dx, dx_2
    cdef double length
    cdef int i, offset
    cdef double * xc
    cdef double * yc

    xc = <double*>malloc(sizeof(double) * N_SAMPLES)
    yc = <double*>malloc(sizeof(double) * N_SAMPLES)

    offset = 0
    for point in range(n_points):
        dx = (bounds[offset + 1] - bounds[offset]) / (N_SAMPLES - 1)
        xc[0] = bounds[offset]
        for i in range(1, N_SAMPLES):
            xc[i] = xc[i - 1] + dx

        c_y_rotated(xc, N_SAMPLES, rw, c, n, a, yc)

        dx_2 = pow(dx, 2)
        length = 0.
        for i in range(1, N_SAMPLES):
            length += sqrt(pow(yc[i] - yc[i-1], 2.) + dx_2)
        lengths[point] = length
        offset += 2

    free(yc)
    free(xc)


@cython.boundscheck(False)
@cython.wraparound(False)
def path_lengths(np.ndarray[np.float_t, ndim=2, mode="c"] bounds not None,
                 double rw, double c, double n, double a,
                 np.ndarray[np.float_t, ndim=1, mode="c"] lengths not None):
    c_path_lengths(&bounds[0, 0], len(bounds), rw, c, n, a, &lengths[0])


cdef void c_nearest_points(double * xy, const int n_points,
                           const double rw, const double c, const double n,
                           const double a, double * nearest) nogil:
    cdef gsl_root_fsolver *s
    cdef gsl_root_fsolver_type *T
    cdef gsl_function F
    cdef double params[5]
    cdef double bounds[2]
    cdef int GSL_CONTINUE = -2
    cdef int status
    cdef double x_of_nearest
    cdef int max_iters = 100
    cdef int iters
    cdef int point
    cdef double x0, y0
    cdef double SIN_A = sin(a)
    cdef double COS_A = cos(a)
    cdef double SIN_MINUS_A = sin(- a)
    cdef double COS_MINUS_A = cos(- a)

    params[2] = rw
    params[3] = c
    params[4] = n

    F.function = &f
    F.params = params

    T = gsl_root_fsolver_brent
    s = gsl_root_fsolver_alloc(T)

    for point in range(0, 2 * n_points, 2):
        x0 = xy[point] * COS_MINUS_A - xy[point + 1] * SIN_MINUS_A
        y0 = xy[point] * SIN_MINUS_A + xy[point + 1] * COS_MINUS_A

        params[0] = x0
        params[1] = y0

        if r(x0, x0, y0, rw, c, n) < 1e-12:
            x_of_nearest = x0
        else:
            find_bounds(params, bounds)

            if nearest_is_lower(params, bounds):
                x_of_nearest = bounds[0]
            else:
                gsl_root_fsolver_set(s, &F, bounds[0], bounds[1])

                status = GSL_CONTINUE
                iters = 0
                while status == GSL_CONTINUE:
                    status = gsl_root_fsolver_iterate(s)
                    x_of_nearest = gsl_root_fsolver_root(s)

                    bounds[0] = gsl_root_fsolver_x_lower(s)
                    bounds[1] = gsl_root_fsolver_x_upper(s)

                    status = gsl_root_test_interval(bounds[0], bounds[1],
                                                    0., 0.001)

                    iters += 1
                    if iters > max_iters:
                        break

        y_of_nearest = centerline_y(x_of_nearest, rw, c, n)

        nearest[point] = x_of_nearest * COS_A - y_of_nearest * SIN_A
        nearest[point + 1] = x_of_nearest * SIN_A + y_of_nearest * COS_A

    gsl_root_fsolver_free(s)


@cython.boundscheck(False)
@cython.wraparound(False)
def nearest_points(np.ndarray[np.float_t, ndim=2, mode="c"] xy not None,
                   double rw, double c, double n, double a,
                   np.ndarray[np.float_t, ndim=2, mode="c"] nearest not None):
    c_nearest_points(&xy[0, 0], len(xy), rw, c, n, a, &nearest[0, 0])

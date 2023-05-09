// ode45.rs

use ndarray::prelude::*;

const C: [f64; 7] = [
    0.0,
    (1.0 / 5.0),
    (3.0 / 10.0),
    (4.0 / 5.0),
    (8.0 / 9.0),
    1.0,
    1.0,
];
const A: [[f64; 7]; 7] = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0],
    [
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
        0.0,
        0.0,
    ],
    [
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0,
        0.0,
        0.0,
    ],
    [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
        0.0,
    ],
];
const BHAT: [f64; 7] = [
    35.0 / 384.0,
    0.0,
    500.0 / 1113.0,
    125.0 / 192.0,
    -2187.0 / 6784.0,
    11.0 / 84.0,
    0.0,
];
const B: [f64; 7] = [
    5179.0 / 57600.0,
    0.0,
    7571.0 / 16695.0,
    393.0 / 640.0,
    -92097.0 / 339200.0,
    187.0 / 2100.0,
    1.0 / 40.0,
];

const MAX_ITER: usize = 1000;

pub fn ode45<F>(
    f: F,
    time_interval: [f64; 2],
    y0: &Array1<f64>,
    tol: f64,
) -> (Vec<f64>, Vec<Array1<f64>>)
where
    F: Fn(f64, &Array1<f64>) -> Array1<f64>,
{
    let mut index = 0usize;

    // Storage for time vector
    let mut times: Vec<f64> = Vec::new();
    times.push(time_interval[0]);

    // Storage for list of values
    let mut values: Vec<Array1<f64>> = Vec::new();
    values.push(y0.clone());

    // Estimate initial timestep
    let dy = f(time_interval[0], y0);
    let d0 = y0.fold(0.0f64, |acc, elem| acc + elem * elem).sqrt();
    let d1 = dy.fold(0.0f64, |acc, elem| acc + elem * elem).sqrt();
    
    // Initial estimate of timestep
    let h0 = if d0 < 1e-5 || d1 < 1e-5 {
        1e-6
    } else { 0.01 * (d0 / d1)};

    // Do euler integration to estimate next time step
    let y1 = y0 + h0 * dy;
    let d2 = y1.fold(0.0f64, |acc, elem| acc + elem * elem).sqrt() / h0;
    let d_max = d2.max(d1);
    let h1 = if d_max <= 1e-15 {
        (1e-6f64).max(h0 * 1e-3)
    } else {
        (0.01 / d_max).powf(1.0 / 5.0)
    };

    // Assign initial timestep based on analysis
    let mut h = h1.min(100.0 * h0);
    println!("initial h: {}", h);

    while times[index] < time_interval[1] {
        for _ in 0..MAX_ITER {
            let mut ks = Vec::with_capacity(7);
            for i in 0..7 {
                let ti = times[index] + h*C[i];
                let yi = values[index].clone()
                    + h * ks.iter()
                        .enumerate()
                        .fold(Array::<f64, _>::zeros(y0.len()), |acc, (j, kj)| {
                            acc + A[i][j] * kj
                        });
                ks.push(f(ti, &yi));
            }

            // Calculate error factor (from Dormand-Prince's paper)
            let e = &ks
                .iter()
                .enumerate()
                .fold(Array::<f64, _>::zeros(y0.len()), |acc, (i, ki)| {
                    acc + (B[i] - BHAT[i]) * ki
                });

            let err = e.fold(0.0f64, |acc, xd| acc + xd * xd).sqrt();

            // Calculate next time step adjustment
            let time_step_factor = ((0.5 * tol * h) / err / (time_interval[1] - time_interval[0])).powf(1.0 / 5.0);
            if 0.9 * time_step_factor < 1.0 {
                h /= 2.0;
            } else if 0.9 * time_step_factor >= 2.0 {
                h *= 2.0;
            } else {
                h *= time_step_factor;
            }

            if err > tol {
                println!("new h: {}", h);
                continue;
            }

            ks.clear();
            for i in 0..7 {
                let ti = times[index] + h*C[i];
                let yi = values[index].clone()
                    + h * ks.iter()
                        .enumerate()
                        .fold(Array::<f64, _>::zeros(y0.len()), |acc, (j, kj)| {
                            acc + A[i][j] * kj
                        });
                ks.push(f(ti, &yi));
            }

            let y_delta = &ks
                .iter()
                .enumerate()
                .fold(Array::<f64, _>::zeros(y0.len()), |acc, (i, ki)| {
                    acc + BHAT[i] * ki
                });

            //println!("y_delta: {}", y_delta);

            values.push(values[index].clone() + h * y_delta);
            times.push(times[index] + h);

            break;
        }
        index += 1;
    }

    (times, values)
}

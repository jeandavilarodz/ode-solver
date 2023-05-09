// rk4.rs

use ndarray::prelude::*;


pub fn rk4<F: Fn(Array1<f64>) -> Array1<f64>>(
    num_steps: usize,
    times: [f64; 2],
    initial_value: Array1<f64>,
    f: F,
) -> (Vec<Array1<f64>>, Array1<f64>) {
    let dt = (times[1] - times[0]) / (num_steps as f64);
    let times = Array::range(0.0, (num_steps + 1) as f64, 1.) * dt + times[0];
    let mut values = Vec::with_capacity(num_steps + 1);
    values.push(initial_value.clone());

    for index in 0..num_steps {
        let k1 = dt * f(values[index].clone());
        let k2 = dt * f(values[index].clone() + (0.5 * k1.clone()));
        let k3 = dt * f(values[index].clone() + (0.5 * k2.clone()));
        let k4 = dt * f(values[index].clone() + k3.clone());
        values.push(values[index].clone() + ((k1 + (2.0 * k2) + (2.0 * k3) + k4) / 6.0));
    }

    (values, times)
}
// 3d_system.rs

use plotters::prelude::*;
use ndarray::prelude::*;
use ode_solver::ode45;

const ALPHA: f64 = 0.10;

// Time series
const YMIN: f64 = -6.0;
const YMAX: f64 = 6.0;

// Phase plot
const XMIN_PHASE: f64 = -6.0;
const XMAX_PHASE: f64 = 6.0;

fn f(_t: f64, xyz: &Array1<f64>) -> Array1<f64> {
    let mut dxdydz = Array::<f64,_>::zeros(xyz.len());
    dxdydz[0] = xyz[1].sin() - ALPHA * xyz[0];
    dxdydz[1] = xyz[2].sin() - ALPHA * xyz[1];
    dxdydz[2] = xyz[0].sin() - ALPHA * xyz[2];
    dxdydz
}

fn main() {
    let y0 = arr1(&[1.0, 2.0, 3.0]);

    let (mut times, mut points) = ode45::ode45(f, [0.0, 1000.0], &y0, 1e-6);
    let points = points.split_off(points.len() / 2);
    let times = times.split_off(times.len() / 2);

    let mut x_values = Vec::with_capacity(points.len());
    let mut y_values = Vec::with_capacity(points.len());
    let mut z_values = Vec::with_capacity(points.len());

    for point in points.iter() {
        x_values.push(point[0]);
        y_values.push(point[1]);
        z_values.push(point[2]);

    }

    //Calculate the mode
    let mut x_values = x_values.split_off(x_values.len() / 2);
    x_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let x1 = x_values[x_values.len() / 2];

    let mut y_values = y_values.split_off(y_values.len() / 2);
    y_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let y1 = y_values[y_values.len() / 2];

    let mut z_values = z_values.split_off(y_values.len() / 2);
    z_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let z1 = z_values[z_values.len() / 2];

    println!("x1: {}, y1: {}, z1: {}", x1, y1, z1);

    let new_initial_conditions = arr1(&[x1, y1, z1]);

    let (mut new_times, mut new_points) = ode45::ode45(f, [0.0, 1000.0], &new_initial_conditions, 1e-6);
    let new_points = new_points.split_off(new_points.len() / 2);
    let new_times = new_times.split_off(new_times.len() / 2);

    // Plot data
    let root = SVGBackend::new("excercise_x.svg", (800, 800)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut plot = ChartBuilder::on(&root)
        .caption(
            format!("x(t) ODE Simulation a = {}", ALPHA),
            ("monospace", (5).percent_height()),
        )
        .y_label_area_size((5).percent_width())
        .x_label_area_size((5).percent_height())
        .margin((1).percent())
        //.build_cartesian_2d((-WIDTH as f32 + X_OFFSET)..(WIDTH as f32 + X_OFFSET) as f32, (-HEIGHT as f32)..(HEIGHT as f32))
        .build_cartesian_2d(500.0..1000.0, YMIN..YMAX)
        .unwrap();
    plot.configure_mesh()
        .draw()
        .unwrap();

    plot.draw_series(
        LineSeries::new(
            points
                .iter()
                .zip(times.iter())
                .map(|(p, &t)| (t, p[0])),
            BLUE,
        )
    )
    .unwrap();

    plot.draw_series(
        LineSeries::new(
            new_points
                .iter()
                .zip(new_times.iter())
                .map(|(p, &t)| (t, p[0])),
            RED,
        )
    )
    .unwrap();

    let root = SVGBackend::new("excercise_y.svg", (800, 800)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut plot = ChartBuilder::on(&root)
        .caption(
            format!("y(t) ODE Simulation a = {}", ALPHA),
            ("monospace", (5).percent_height()),
        )
        .y_label_area_size((5).percent_width())
        .x_label_area_size((5).percent_height())
        .margin((1).percent())
        //.build_cartesian_2d((-WIDTH as f32 + X_OFFSET)..(WIDTH as f32 + X_OFFSET) as f32, (-HEIGHT as f32)..(HEIGHT as f32))
        .build_cartesian_2d(500.0..1000.0, YMIN..YMAX)
        .unwrap();
    plot.configure_mesh()
        .draw()
        .unwrap();

    plot.draw_series(
        LineSeries::new(
            points
                .iter()
                .zip(times.iter())
                .map(|(p, &t)| (t, p[1])),
            BLUE,
        )
    )
    .unwrap();

    plot.draw_series(
        LineSeries::new(
            new_points
                .iter()
                .zip(new_times.iter())
                .map(|(p, &t)| (t, p[1])),
            RED,
        )
    )
    .unwrap();

    let root = SVGBackend::new("excercise_phase.svg", (800, 800)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut plot = ChartBuilder::on(&root)
        .caption(
            format!("XZ-Phase ODE Simulation a = {}", ALPHA),
            ("monospace", (5).percent_height()),
        )
        .y_label_area_size((5).percent_width())
        .x_label_area_size((5).percent_height())
        .margin((1).percent())
        //.build_cartesian_2d((-WIDTH as f32 + X_OFFSET)..(WIDTH as f32 + X_OFFSET) as f32, (-HEIGHT as f32)..(HEIGHT as f32))
        .build_cartesian_2d(XMIN_PHASE..XMAX_PHASE, XMIN_PHASE..XMAX_PHASE)
        .unwrap();
    plot.configure_mesh()
        .draw()
        .unwrap();

    plot.draw_series(
        LineSeries::new(
            points
                .iter()
                .map(|p| (p[0], p[2])),
            BLUE,
        )
    )
    .unwrap();

    plot.draw_series(
        LineSeries::new(
            new_points
                .iter()
                .map(|p| (p[0], p[2])),
            RED,
        )
    )
    .unwrap();

}

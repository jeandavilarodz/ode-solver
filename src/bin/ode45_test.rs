// main.rs

use plotters::prelude::*;
use ndarray::prelude::*;
use ode_solver::ode45;

fn f(t: f64, y: &Array1<f64>) -> Array1<f64> {
    let mut ret = Array::<f64,_>::zeros(y.len());
    ret[0] = y[0] * (2.0 - t) * t + t - 1.0;
    ret
}

fn main() {
    let y0 = arr1(&[1.0]);

    let (times, points) = ode45::ode45(f, [0.0, 5.0], &y0, 0.1, 1e-6);

    println!("{:?}", times);

    // Plot data
    let root = SVGBackend::new("ode45.svg", (1280, 720)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut plot = ChartBuilder::on(&root)
        .caption(
            "ODE Simulation",
            ("monospace", (5).percent_height()),
        )
        .y_label_area_size((5).percent_width())
        .x_label_area_size((5).percent_height())
        .margin((1).percent())
        //.build_cartesian_2d((-WIDTH as f32 + X_OFFSET)..(WIDTH as f32 + X_OFFSET) as f32, (-HEIGHT as f32)..(HEIGHT as f32))
        .build_cartesian_2d(0.0..5.0, 0.0..3.0)
        .unwrap();
    plot.configure_mesh()
        .draw()
        .unwrap();

    plot.draw_series(
        LineSeries::new(
            points
                .into_iter()
                .zip(times.iter())
                .map(|(p, t)| (*t, p[0])),
            BLUE.stroke_width(2),
        ).point_size(2)
    )
    .unwrap();

}

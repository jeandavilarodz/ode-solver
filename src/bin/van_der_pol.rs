// main.rs

use plotters::prelude::*;
use ndarray::prelude::*;
use ode_solver::ode45;

fn f(_t: f64, p: &Array1<f64>) -> Array1<f64> {
    let mut dxdy = Array::<f64,_>::zeros(p.len());
    //ret[0] = y[0] * (2.0 - t) * t + t - 1.0;
    dxdy[0] = p[1];
    dxdy[1] = 2.0 * (1.0 - p[0] * p[0]) * p[1] - p[0];
    dxdy
}

fn main() {
    let y0 = arr1(&[0.5, 0.0]);

    let (times, points) = ode45::ode45(f, [0.0, 20.0], &y0, 1e-6);

    // Plot data
    let root = SVGBackend::new("van_der_pol.svg", (800, 800)).into_drawing_area();
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
        .build_cartesian_2d(0.0..20.0, -3.0..3.0)
        .unwrap();
    plot.configure_mesh()
        .draw()
        .unwrap();

    plot.draw_series(
        LineSeries::new(
            points
                .iter()
                .zip(times.iter())
                .map(|(p, t)| (*t, p[0])),
            BLUE,
        ).point_size(1)
    )
    .unwrap();

    let root = SVGBackend::new("van_der_pol_phase.svg", (800, 800)).into_drawing_area();
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
        .build_cartesian_2d(-4.0..4.0, -4.0..4.0)
        .unwrap();
    plot.configure_mesh()
        .draw()
        .unwrap();
    plot.draw_series(
        LineSeries::new(
            points
                .iter()
                .map(|p| (p[0], p[1])),
            BLUE,
        ).point_size(1)
    )
    .unwrap();

}

// main.rs

use plotters::prelude::*;
use ndarray::prelude::*;
mod rk4;


fn main() {
    let f = |x: Array1<f64>| {
        let mut dxdt = Array::zeros(x.len());
        dxdt[0] = 2.0 * x[0];
        dxdt
    };

    let x: Array1<f64> = Array1::from_vec(vec![1.0 / 4.0,]);

    let (points, t) = rk4::rk4(64, [0.0, 2.0], x, &f);

    println!("{:?}", points.last());

    // Plot data
    let root = SVGBackend::new("ode-solution-map.svg", (1280, 720)).into_drawing_area();
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
        .build_cartesian_2d(0.0..2.0, 0.0..14.0)
        .unwrap();
    plot.configure_mesh()
        .draw()
        .unwrap();

    plot.draw_series(
        LineSeries::new(
            points
                .into_iter()
                .zip(t.iter())
                .map(|(p, t)| (*t, p[0])),
            BLUE.stroke_width(2),
        ).point_size(2)
    )
    .unwrap();

    plot.draw_series(
        LineSeries::new(
            t.iter().map(|t| (*t, 0.25 * ( 2.0 * t).exp())),
            GREEN.stroke_width(1),
        )
    )
    .unwrap();
}

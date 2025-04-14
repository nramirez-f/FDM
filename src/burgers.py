import plotly.graph_objects as go
import numpy as np

def plot_method_of_characteristics(x0:float, xf:float, nx:int, T:float, nt:int, f):
    """
    Plots the solution of the Burgers' equation using the method of characteristics.

    Parameters:
    -----------
    x0 : float
        The initial spatial coordinate (left boundary of the domain).
    
    xf : float
        The final spatial coordinate (right boundary of the domain).
    
    nx : int
        The number of spatial grid points.
    
    T : float
        The total simulation time.
    
    nt : int
        The number of time steps.
    
    f : function
        The initial condition function, which defines the initial profile of the solution.

    Returns:
    --------
    None
        The function generates and displays a plot but does not return any values.
    """
    x = np.linspace(x0, xf, nx)
    xi = np.copy(x)
    dt = T / nt

    # Figura inicial
    fig = go.Figure(
        data=[
            go.Scatter(
                x=x,
                y=f(x),
                mode="lines",
                showlegend=False,
            )
        ]
    )

    iteration_frames = []

    for k in range(nt + 1):
        t = k * dt
        x = xi + f(xi) * t
        iteration_frames.append(
            go.Frame(
                data=[
                    go.Scatter(
                        x=x,
                        y=f(xi),
                        mode="lines",
                        showlegend=False,
                    )
                ],
                name=str(t),
            )
        )

    fig.update(frames=iteration_frames)

    fig.update_layout(
        title="Burgers simulation by Method of Characteristics",
        yaxis=dict(
            range=[f(xi).min() - 0.1, f(xi).max() + 0.1],
            title="u",
        ),
        xaxis=dict(
            range=[min(x[0], xi[0]), max(x[-1], xi[-1])],
            title="x",
        ),
    )

    # Controles de animación
    fig.update_layout(
        updatemenus=[
            {
                "buttons": [
                    {
                        "args": [None, {"frame": {"duration": 200, "redraw": True}, "fromcurrent": True}],
                        "label": "▶ Play",
                        "method": "animate",
                    },
                    {
                        "args": [[None], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate", "transition": {"duration": 0}}],
                        "label": "❚❚ Pause",
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 87},
                "showactive": False,
                "type": "buttons",
                "x": 0.1,
                "xanchor": "right",
                "y": 0,
                "yanchor": "top",
            }
        ]
    )

    fig.update_layout(
        sliders=[
            {
                "steps": [
                    {
                        "args": [[str(iterValue)], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}],
                        "label": str(iterValue),
                        "method": "animate",
                    }
                    for iterValue in np.arange(0, T + dt, dt)
                ],
                "x": 0.1,
                "len": 0.9,
                "xanchor": "left",
                "y": -0.2,
                "yanchor": "top",
            }
        ]
    )

    fig.show()
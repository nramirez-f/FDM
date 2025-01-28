# -*- coding: utf-8 -*-
from netCDF4 import Dataset
import plotly.graph_objects as go


class ncv:
    def __init__(self, filepath):
        """
        Initialize the view with a given file.
        
        Parameters:
        - filepath: Path to the NetCDF file.
        """
        # Abrir el archivo NetCDF
        rootgrp = Dataset(filepath, "r")
        
        # Extraer las variables del archivo
        times = rootgrp.variables["time"][:]  # Tiempo
        x = rootgrp.variables["x"][:]         # Coordenada espacial
        u = rootgrp.variables["u"][:]         # Variable u (simulación)

        # Calcular el valor mínimo y máximo de u a través de todos los tiempos y espacios
        u_min = u.min()
        u_max = u.max()

        # Agregar un pequeño margen al rango
        margin = 0.05  # Margen de 5% por encima y por debajo del rango
        u_range_min = u_min - margin * (u_max - u_min)
        u_range_max = u_max + margin * (u_max - u_min)
        
        # Crear la figura
        fig = go.Figure(
            data=[go.Scatter(x=x, y=u[0, :], mode="lines", name=f"t = {times[0]:.2f} s")],
            layout=go.Layout(
                title="Evolución de la variable u a lo largo del tiempo",
                xaxis_title="Espacio (x)",
                yaxis_title="Valor de u",
                showlegend=True,
                updatemenus=[
                    dict(
                        type="buttons",
                        x=0.1,
                        xanchor="right",
                        y=1.15,
                        yanchor="top",
                        buttons=[
                            dict(
                                label="Play",
                                method="animate",
                                args=[None, dict(frame=dict(duration=100, redraw=True), fromcurrent=True)]
                            ),
                            dict(
                                label="Pause",
                                method="animate",
                                args=[[None], dict(frame=dict(duration=0, redraw=True), mode="immediate", transition=dict(duration=0))]
                            )
                        ]
                    )
                ],
                sliders=[dict(
                    yanchor="top",
                    xanchor="left",
                    currentvalue=dict(
                        visible=True,
                        prefix="Tiempo: ",
                        font=dict(size=20),
                        offset=20
                    ),
                    steps=[]
                )],
                yaxis=dict(
                    range=[u_range_min, u_range_max]  # Establecer el rango del eje y
                )
            )
        )

        # Crear las frames para la animación
        frames = []
        for t_index in range(1, len(times)):
            frame = go.Frame(
                data=[go.Scatter(x=x, y=u[t_index, :], mode="lines", name=f"t = {times[t_index]:.2f} s")],
                name=f"t = {times[t_index]:.2f} s"
            )
            frames.append(frame)
        
        # Agregar las frames a la figura
        fig.frames = frames

        # Crear los pasos para el slider
        steps = []
        for t_index in range(len(times)):
            step = dict(
                method="animate",
                args=[
                    [f"t = {times[t_index]:.2f} s"],
                    dict(
                        frame=dict(duration=100, redraw=True),
                        mode="immediate",
                        transition=dict(duration=0)
                    )
                ],
                label=f"{times[t_index]:.2f} s"
            )
            steps.append(step)
        
        # Actualizar la configuración del slider
        fig.layout.sliders[0]["steps"] = steps

        # Mostrar la visualización
        fig.show()

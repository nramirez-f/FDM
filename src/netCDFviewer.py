# -*- coding: utf-8 -*-
from netCDF4 import Dataset
import plotly.graph_objects as go




class NCV:
    """
    Object of netCDF to visualize it.
    """

    rootgrp = None

    def __init__(self, filepath):
        """
        Initialize the netCDF object with the path of the file.
    
        - filepath: Path to the NetCDF file.
        """
    
        if (filepath is None or filepath ==""):
            raise Exception("filepath must be distinct of empty string.")
        else:
            self.rootgrp = Dataset(filepath, "r")


    def shape(self, iterPos, iterName = "time", varName = "u", dimName = "x"):
        """
            Static plot at one iterable position.

            - itername: Name of the iterable dimension
            - 
        """

        iterDomain = self.rootgrp.variables[iterName][:]

        if iterPos in iterDomain:
            fig = go.Figure()

            fig.add_trace(go.Scatter(x=self.rootgrp.variables[dimName][:], y=self.rootgrp.variables[varName][:, iterPos],
                                    mode='lines',
                                    name=varName+"("+ dimName + ", " + str(iterPos) +")",
                                    showlegend=True))
            
            fig.show()
        else:
            raise Exception("iterPos not in its posible values.")

    def playShape(self, iterInit, iterName="time", varName="u", dimName="x"):
        """
            Animated plot from iterInit to the end of iterDomain.

            - iterInit: Initial value for the iteration
            - iterName: Name of the iterable dimension
            - varName: Name of the variable
            - dimName: Name of the dimension
        """
        fullIterDomain = self.rootgrp.variables[iterName][:]

        if iterInit not in fullIterDomain:
            raise Exception("iterInit not in its possible values.")
        else:
            iterDomain = [i for i in fullIterDomain if i >= iterInit]

            fig = go.Figure(
                data=[
                    go.Scatter(
                        x=self.rootgrp.variables[dimName][:],
                        y=self.rootgrp.variables[varName][:][:, iterDomain.index(iterInit)],
                        mode="lines",
                        name=f"{varName}({dimName}, {iterInit})"
                    )
                ]
            )

            eps_y = 0.1
            fig.update_layout(
                yaxis=dict(range=[self.rootgrp.variables[varName][:].min() - eps_y, self.rootgrp.variables[varName][:].max() + eps_y])
            )

            frames = [
                go.Frame(
                    data=[
                        go.Scatter(
                            x=self.rootgrp.variables[dimName][:],
                            y=self.rootgrp.variables[varName][:][:, iterDomain.index(iterValue)],
                            mode="lines",
                            name=f"{varName}({dimName}, {iterValue})"
                        )
                    ],
                    name=str(iterValue)
                )
                for iterValue in iterDomain
            ]

            fig.update(frames=frames)

            fig.update_layout(
                updatemenus=[
                    {
                        "buttons": [
                            {
                                "args": [None, {"frame": {"duration": 200, "redraw": True}, "fromcurrent": True}],
                                "label": "▶ Play",
                                "method": "animate"
                            },
                            {
                                "args": [[None], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate", "transition": {"duration": 0}}],
                                "label": "❚❚ Pause",
                                "method": "animate"
                            }
                        ],
                        "direction": "left",
                        "pad": {"r": 10, "t": 87},
                        "showactive": False,
                        "type": "buttons",
                        "x": 0.1,
                        "xanchor": "right",
                        "y": 0,
                        "yanchor": "top"
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
                                "method": "animate"
                            }
                            for iterValue in iterDomain
                        ],
                        "x": 0.1,
                        "len": 0.9,
                        "xanchor": "left",
                        "y": -0.2,
                        "yanchor": "top"
                    }
                ]
            )

            fig.show()


            
    def close(self):
        self.rootgrp.close()

    ### WITH CONTEXT ###

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

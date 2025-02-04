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
    
    ### BASIC ###
    def scatter(self, dimName, varName, iterValue, legend=""):
        sct = go.Scatter(
                x=self.rootgrp.variables[dimName][:],
                y=self.rootgrp.variables[varName][:][:, iterValue],
                mode="lines",
                showlegend=bool(legend),
                name=legend,
            )

        return sct
    
    ### PLOT ###
    def shape(self, iterValue, iterName = "time", varName = "u", dimName = "x"):
        """
            Static plot at one iterable position.

            - itername: Name of the iterable dimension
            - 
        """

        iterDomain = self.rootgrp.variables[iterName][:]

        if iterValue in iterDomain:
            fig = go.Figure()

            fig.add_trace(self.scatter(dimName, varName, iterValue, f"{varName}({dimName}, {iterValue})"))
            
            # Layout
            fig.update_layout(
                    title=dict(
                        text=self.rootgrp.description
                    ),
                    xaxis=dict(
                        title=dict(
                            text=dimName
                        )
                    ),
                    yaxis=dict(
                        title=dict(
                            text=varName
                        )
                    ),
            )

            fig.show()
        else:
            raise Exception("iterValue not in its posible values.")

    def playShape(self, iterInit, delay=100, iterName="time", varName="u", dimName="x"):
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
                    self.scatter(dimName, varName, iterDomain.index(iterInit))
                ]
            )

            # Frames
            fig.update(
                frames=[
                    go.Frame(
                        data=[
                            self.scatter(dimName, varName, iterDomain.index(iterValue))
                        ],
                        name=str(iterValue)
                    )
                    for iterValue in iterDomain
                ]
            )

            # Layout
            fig.update_layout(
                title=dict(text=self.rootgrp.description),
            )

            eps_y = 0.1
            fig.update_layout(
                yaxis=dict(range=[self.rootgrp.variables[varName][:].min() - eps_y, self.rootgrp.variables[varName][:].max() + eps_y],
                           title=dict(text=varName)),
            )
            
            fig.update_layout(
                    xaxis=dict(
                        range = [self.rootgrp.variables[dimName][0], self.rootgrp.variables[dimName][-1]],
                        title=dict(text=dimName),
                    ),
            )

            # Animation Controlsb
            fig.update_layout(
                updatemenus=[
                    {
                        "buttons": [
                            {
                                "args": [None, {"frame": {"duration": delay, "redraw": True}, "fromcurrent": True}],
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

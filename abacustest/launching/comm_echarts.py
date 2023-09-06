import numpy as np
def trans_nan_to_none(x):
    if x is None:
        return None
    if isinstance(x, float) and np.isnan(x):
        return None
    return x


def get_bar_option(title, x, y, x_type="category", y_type="value"):
    y_used = [trans_nan_to_none(i) for i in y]
    if y_type == "log":
        # need set 0 to None and set negative to abs
        for ix, v in enumerate(y_used):
            if v in [0, None]:
                y_used[ix] = None
            elif v < 0:
                y_used[ix] = -v

    option = {
        "title": {
            "text": title,
            "textStyle": {
                "color": '#333',
                "fontSize": 18,
                "fontWeight": 'bold'
            }
        },
        "toolbox": {
            "show": True,
            "feature": {
                "dataZoom": {},
                "dataView": {"readOnly": False},
                "magicType": {"type": ['line', 'bar']},
                "restore": {},
                "saveAsImage": {}
            }
        },
        "tooltip": {
            "trigger": 'axis',
            "axisPointer": {
                "type": 'shadow'
            }
        },
        "grid": {
            "left": '3%',
            "right": '4%',
            "bottom": '3%',
            "containLabel": True
        },
        "xAxis": [
            {
                "type": x_type,
                "data": x,
                "axisTick": {
                    "alignWithLabel": True
                },
            }
        ],
        "yAxis": [
            {
                "type": y_type,
                "scale": True
            }
        ],
        "series": [
            {
                "name": 'Direct',
                "type": 'line',
                "barWidth": '60%',
                "data": y_used
            }
        ]
    }
    return option


def produce_multiple_y(title, x, y_list, legend_list, x_type="category", y_type="value"):
    series = []
    y_list = [[trans_nan_to_none(i) for i in y] for y in y_list]
    for iy, y in enumerate(y_list):
        series.append({
            "name": legend_list[iy],
            "type": 'line',
            "barWidth": '60%',
            "data": y.copy(),
            "gridIndex": 0,
            "yAxisIndex": 0,
            "xAxisIndex": 0,
        })
    grid = [{
            "left": '3%',
            "right": '4%',
            "bottom": '3%',
            "top": "25%",
            "containLabel": True
            }]
    xAxis= [
            {
                "type": x_type,
                "data": x,
                "axisTick": {
                    "alignWithLabel": True
                },
                "gridIndex": 0,
            }
        ]
    yAxis = [
            {
                "type": "value",
                "scale": True,
                "gridIndex": 0,
            }
        ]
    if y_type == "log":
        grid = [{
            "left": '3%',
            "right": '55%',
            "bottom": '3%',
            "top": "25%",
            "containLabel": True
        }, {
            "left": '50%',
            "right": '3%',
            "bottom": '3%',
            "top": "25%",
            "containLabel": True
        }]
        xAxis.append({
                "type": x_type,
                "data": x,
                "axisTick": {
                    "alignWithLabel": True
                },
                "gridIndex": 1,
            })
        yAxis.append({
                "type": "log",
                "scale": True,
                "gridIndex": 1,
            })

        # need set 0 to None and set negative to abs
        log_series = []
        for iy, y in enumerate(y_list):
            log_series.append({
                "name": legend_list[iy],
                "type": 'line',
                "barWidth": '60%',
                "data": y.copy(),
                "gridIndex": 1,
                "yAxisIndex": 1,
                "xAxisIndex": 1,
            })
            for ix, v in enumerate(y):
                if v in [0, None]:
                    log_series[iy]["data"][ix] = None
                elif v < 0:
                    log_series[iy]["data"][ix] = -v
        series += log_series

    option = {
        "title": {
            "text": title,
            "textStyle": {
                "color": '#333',
                "fontSize": 18,
                "fontWeight": 'bold'
            }
        },
        "toolbox": {
            "show": True,
            "feature": {
                "dataZoom": {},
                "dataView": {"readOnly": False},
                "magicType": {"type": ['line', 'bar']},
                "restore": {},
                "saveAsImage": {}
            }
        },
        "tooltip": {
            "trigger": 'axis',
            "axisPointer": {
                "type": 'shadow'
            }
        },
        "grid": grid,
        "xAxis": xAxis,
        "yAxis": yAxis,
        "series": series,
        "legend": {"data": legend_list,
                   "top": "10%"}
    }
    return option

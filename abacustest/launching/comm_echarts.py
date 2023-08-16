
def get_bar_option(title,x,y,x_type="category",y_type="value"):
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
                "data": y
            }
        ]
    }
    return option

def produce_multiple_y(title,x,y_list,legend_list,x_type="category",y_type="value"):
    series = []
    for iy,y in enumerate(y_list):
        series.append({
                "name": legend_list[iy],
                "type": 'line',
                "barWidth": '60%',
                "data": y
            })
        
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
        "series": series,
        "legend": {"data": legend_list}
    }
    return option
    

def get_bar_option(title,x,y):
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
                "type": 'category',
                "data": x,
                "axisTick": {
                    "alignWithLabel": True
                },
            }
        ],
        "yAxis": [
            {
                "type": 'value',
                "scale": True
            }
        ],
        "series": [
            {
                "name": 'Direct',
                "type": 'bar',
                "barWidth": '60%',
                "data": y
            }
        ]
    }
    return option

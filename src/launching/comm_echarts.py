
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
                "type": 'value'
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

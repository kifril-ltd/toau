import json
from json import JSONEncoder

from django.http import JsonResponse, FileResponse, HttpResponse
from django.shortcuts import render
from django.views import View
import numpy

from stp.utils.StpOptimizer import StpOptimizer


class MainView(View):
    def get(self, request):
        return render(request, 'main.html')

class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)

def calculateOpimization(request):
    inputData = request.POST.get('data')
    dataDict = json.loads(inputData)

    stpOptimizer = StpOptimizer(dataDict)

    sol = stpOptimizer.optimizeStp()
    sol_non = stpOptimizer.optimizeNonStp()

    rif = stpOptimizer.calcRIF()
    ricr = stpOptimizer.calcRICR()
    stat = stpOptimizer.pValueStat()

    context = {
        'stp_sol': sol,
        'non_stp_sol': sol_non,
        'rif': rif,
        'ricr': ricr,
        'p_solutions': stat,
        'products_number': dataDict['products_number'],
        'resources_number': dataDict['resources_number'],
    }

    context = json.dumps(context, cls=NumpyArrayEncoder)

    return HttpResponse(context)

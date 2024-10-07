from django.shortcuts import render
# Create your views here.
from django.http import HttpResponse
from sympy import Symbol, sympify
import pandas as pd

def index(request):
    return render(request, 'parameters.html')

def result(request):
    method = str(request.GET.get('method'))
    functionStr = str(request.GET.get('functionStr'))
    yInit = str(request.GET.get('yInit'))
    yInit = sympify(yInit)
    hStep = str(request.GET.get('hStep'))
    hStep = sympify(hStep)
    tVal = str(request.GET.get('tVal'))
    tVal = sympify(tVal)
    t = Symbol('t')
    y = Symbol('y')
    functionExpr = sympify(functionStr)
    if method == "Euler":
        approximation = euler(functionExpr, yInit, hStep,0,tVal)[0]
        processH = euler(functionExpr, yInit, hStep,0,tVal)[1]
    elif method == "Improved Euler":
        approximation = improvedEuler(functionExpr, yInit, hStep,0,tVal)[0]
        processH = improvedEuler(functionExpr, yInit, hStep,0,tVal)[1]
    elif method == "Runge-Kutta":
        approximation = rungeKutta(functionExpr, yInit, hStep,0,tVal)[0]
        processH = rungeKutta(functionExpr, yInit, hStep,0,tVal)[1]
    return render(request, 'result.html', {'method': method, 'approximation':approximation, 'processH':processH})

def euler(expression, currY, step, tInit, tFin):
    process = pd.DataFrame(data = {'n':[0],'t':[tInit],'y':[currY]})
    count = 0
    while (tInit < tFin):
        t = Symbol('t')
        y = Symbol('y')
        h= Symbol('h')
        euler_expression = y+h*(expression)
        currY = euler_expression.subs([(y,currY),(t, tInit),(h,step)])
        tInit += step
        count+=1
        process = process._append({'n': int(count), 't': tInit, 'y': float(currY)},ignore_index=True)
    processH = process.to_html
    return [currY, processH]

def improvedEuler(expression, currY, step, tInit, tFin):
    process = pd.DataFrame(data = {'n':[0],'t':[tInit],'y':[currY]})
    count = 0
    while (tInit < tFin):
        t = Symbol('t')
        y = Symbol('y')
        h= Symbol('h')
        expression2 = expression.subs([(t,t+h),(y,y+(h*(expression)))])
        impEulerExpression = y+((h/2)*(expression+expression2))
        currY = impEulerExpression.subs([(y,currY),(t, tInit),(h,step)])
        tInit += step
        count+=1
        process = process._append({'n': int(count), 't': tInit, 'y': float(currY)},ignore_index=True)   
    processH = process.to_html
    return [currY,processH]

def rungeKutta(expression, currY, step, tInit, tFin):
    process = pd.DataFrame(data = {'n':[0],'t':[tInit],'y':[currY]})
    count = 0
    while (tInit < tFin):
        t = Symbol('t')
        y = Symbol('y')
        h= Symbol('h')
        expression2 = expression.subs([(t,t+(.5*h)),(y,y+(.5*h)*(expression))])
        expression3 = expression.subs([(t,t+(.5*h)),(y,y+(.5*h)*(expression2))])
        expression4 = expression.subs([(t,t+h),(y,y+h*(expression3))])
        rungeKuttaExpression = y+((h/6)*(expression+(2*expression2)+(2*expression3)+expression4))
        currY = rungeKuttaExpression.subs([(y,currY),(t, tInit),(h,step)])
        tInit += step
        count+=1
        process = process._append({'n': int(count), 't': tInit, 'y': float(currY)},ignore_index=True)   
    processH = process.to_html
    return [currY, processH]

t = Symbol('t')
y = Symbol('y')
print(improvedEuler(sympify(2+t-y),1,.025,0,.4)[1])
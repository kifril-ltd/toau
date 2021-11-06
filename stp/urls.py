from django.urls import path

from . import views

urlpatterns = [
    path('', views.MainView.as_view()),
    path('calc', views.calculateOpimization, name='calc'),
]
/*
 * 文件名: fittingchart.h
 * 文件作用: 拟合绘图管理类头文件
 * 修改记录:
 * 1. [优化] 增加 initializeCharts 对坐标轴的预配置。
 * 2. [优化] plotAll 增加 autoScale 参数，默认为 false。
 */

#ifndef FITTINGCHART_H
#define FITTINGCHART_H

#include <QObject>
#include <QVector>
#include "mousezoom.h"
#include "fittingdatadialog.h"

class FittingChart : public QObject
{
    Q_OBJECT
public:
    explicit FittingChart(QObject *parent = nullptr);

    // 初始化图表指针，并进行一次性的坐标轴风格配置
    void initializeCharts(MouseZoom* logLog, MouseZoom* semiLog, MouseZoom* cartesian);

    void setObservedData(const QVector<double>& t, const QVector<double>& deltaP,
                         const QVector<double>& deriv, const QVector<double>& rawP);

    void setSettings(const FittingDataSettings& settings);

    // 绘制所有曲线
    // autoScale: 是否自动缩放视图。默认为 false (保持用户当前视图)
    void plotAll(const QVector<double>& t_model, const QVector<double>& p_model, const QVector<double>& d_model, bool isModelValid, bool autoScale = false);

    // 绘制抽样点
    void plotSampledPoints(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d);

    double getCalculatedInitialPressure() const;

private:
    MouseZoom* m_plotLogLog;
    MouseZoom* m_plotSemiLog;
    MouseZoom* m_plotCartesian;

    QVector<double> m_obsT;
    QVector<double> m_obsDeltaP;
    QVector<double> m_obsDeriv;
    QVector<double> m_obsRawP;

    FittingDataSettings m_settings;
    double m_calculatedPi;

    // 内部绘图函数
    void plotLogLog(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel, bool autoScale);
    void plotSemiLog(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel, bool autoScale);
    void plotCartesian(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel, bool autoScale);

    double calculateHornerPressure();
    void showResultOnLogPlot();
};

#endif // FITTINGCHART_H

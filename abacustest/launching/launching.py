import traceback
from dp.launching.cli import SubParser,run_sp_and_exit,default_exception_handler
from abacustest.launching import (model_normal,
 model_selfDefine,
 model_summary,
 model_postdft,
 model_predft,
 model_advanced,
 model_expert,
 model_report,
 model_reuse,
 model_phonon,
 model_fdforce,
 model_fdstress
 )


def to_parser():
    return {
       # "Normal":SubParser(model_normal.NormalModel,model_normal.NormalModelRunner,"Simple rundft and postdft"),
       # "Advanced":SubParser(model_advanced.AdvancedModel,model_advanced.AdvancedModelRunner,"Complete rundft and postdft"),
       # "Expert":SubParser(model_expert.ExpertModel,model_expert.ExpertModelRunner,"Complete predft, rundft and postdft"),
       # "Predft":SubParser(model_predft.PredftModel,model_predft.PredftModelRunner,"Only Predft"),
       # "Postdft":SubParser(model_postdft.PostdftModel,model_postdft.PostdftModelRunner,"Only Postdft"),
       # "SettingFile":SubParser(model_selfDefine.SelfDefineModel,model_selfDefine.SelfDefineModelRunner,"run self-defined model"),
        #"NormalDatasets":SubParser(model_normal.NormalDatasetsModel,model_normal.NormalModelRunner,"Simple rundft and postdft (use launching datasets as input)"),
        #"AdvancedDatasets":SubParser(model_advanced.AdvancedDatasetsModel,model_advanced.AdvancedModelRunner,"Complete rundft and postdft (use launching datasets as input)"),
       # "ExpertDatasets":SubParser(model_expert.ExpertDatasetsModel,model_expert.ExpertModelRunner,"Complete predft, rundft and postdft (use launching datasets as input)"),
        #"PredftDatasets":SubParser(model_predft.PredftDatasetsModel,model_predft.PredftModelRunner,"Only Predft (use launching datasets as input)"),
        #"PostdftDatasets":SubParser(model_postdft.PostdftDatasetsModel,model_postdft.PostdftModelRunner,"Only Postdft (use launching datasets as input)"),
       # "SettingFileDatasets":SubParser(model_selfDefine.SelfDefineDatasetsModel,model_selfDefine.SelfDefineModelRunner,"run self-defined model(use launching datasets as input)"),
        #"UploadDatasets":SubParser(model_uploadDataset.UplaodDatasetModel,model_uploadDataset.UplaodDatasetModelRunner,"upload datasets to datahub"),
       # "Report":SubParser(model_report.ReportModel,model_report.ReportModelRunner,"report metrics.json"),
        "00-Reuse":SubParser(model_reuse.ReuseModel,model_reuse.ReuseModelRunner,"reuse model, use model to run other datasets"),
        "01-Report":SubParser(model_report.ReportModel,model_report.ReportModelRunner,"Show abacustest.html, supermetrics and metrics"),       
        "02-Phonon": SubParser(model_phonon.PhononModel,model_phonon.PhononModelRunner,"Calculate phonon"),    
        "03-FDForce":  SubParser(model_fdforce.FDForceModel,model_fdforce.FDForceModelRunner,"Do finite difference force calculation. Need info.txt file in each example inputs."), 
        "04-FDStress":  SubParser(model_fdstress.FDStressModel,model_fdstress.FDStressModelRunner,"Do finite difference stress calculation."),  
        
        #"Summary":SubParser(model_summary.SummaryModel,model_summary.SummaryModelRunner,"summary abacustest results"),

    }

if __name__ == '__main__':
    run_sp_and_exit(to_parser(), 
                    description="abacustest", 
                    version="0.1.0",
                    exception_handler=default_exception_handler)
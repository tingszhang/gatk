-- MySQL dump 10.13  Distrib 5.7.19, for osx10.12 (x86_64)
--
-- Host: mysql-prd2.broadinstitute.org    Database: DSPRegressionTesting
-- ------------------------------------------------------
-- Server version	5.7.17-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `AnalysisInfo`
--

DROP TABLE IF EXISTS `AnalysisInfo`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `AnalysisInfo` (
  `idAnalysisInfo` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `version` varchar(1024) COLLATE utf8mb4_unicode_ci NOT NULL,
  `name` varchar(512) COLLATE utf8mb4_unicode_ci NOT NULL,
  `wdl` varchar(1024) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `wdlChecksum` char(32) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`idAnalysisInfo`),
  UNIQUE KEY `idAnalysisInfo_UNIQUE` (`idAnalysisInfo`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `AnalysisInfo`
--

LOCK TABLES `AnalysisInfo` WRITE;
/*!40000 ALTER TABLE `AnalysisInfo` DISABLE KEYS */;
/*!40000 ALTER TABLE `AnalysisInfo` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `AnalysisOutputFiles`
--

DROP TABLE IF EXISTS `AnalysisOutputFiles`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `AnalysisOutputFiles` (
  `analysis` int(11) unsigned NOT NULL,
  `outputFile` int(11) unsigned NOT NULL,
  PRIMARY KEY (`analysis`,`outputFile`),
  KEY `AnalysisOutputFiles_outputFileID_idx` (`outputFile`),
  CONSTRAINT `AnalysisOutputFiles_analysisID` FOREIGN KEY (`analysis`) REFERENCES `AnalysisRuns` (`idAnalysisRun`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `AnalysisOutputFiles_outputFileID` FOREIGN KEY (`outputFile`) REFERENCES `OutputFiles` (`idOutputFiles`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `AnalysisOutputFiles`
--

LOCK TABLES `AnalysisOutputFiles` WRITE;
/*!40000 ALTER TABLE `AnalysisOutputFiles` DISABLE KEYS */;
/*!40000 ALTER TABLE `AnalysisOutputFiles` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `AnalysisRuns`
--

DROP TABLE IF EXISTS `AnalysisRuns`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `AnalysisRuns` (
  `idAnalysisRun` int(11) unsigned NOT NULL,
  `analysisInfoID` int(11) unsigned NOT NULL,
  `scenarioOutputForComparison` int(11) unsigned DEFAULT NULL,
  `configuration` json DEFAULT NULL,
  PRIMARY KEY (`idAnalysisRun`),
  UNIQUE KEY `idAnalysis_UNIQUE` (`idAnalysisRun`),
  KEY `Analysis_scenarioOutputForComparison_idx` (`scenarioOutputForComparison`),
  KEY `AnalysisRuns_analysisInfoID_idx` (`analysisInfoID`),
  CONSTRAINT `AnalysisRuns_analysisInfoID` FOREIGN KEY (`analysisInfoID`) REFERENCES `AnalysisInfo` (`idAnalysisInfo`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `AnalysisRuns_scenarioOutputForComparison` FOREIGN KEY (`scenarioOutputForComparison`) REFERENCES `OutputFiles` (`idOutputFiles`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `AnalysisRuns`
--

LOCK TABLES `AnalysisRuns` WRITE;
/*!40000 ALTER TABLE `AnalysisRuns` DISABLE KEYS */;
/*!40000 ALTER TABLE `AnalysisRuns` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `FileTypes`
--

DROP TABLE IF EXISTS `FileTypes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `FileTypes` (
  `idFileTypes` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `type` varchar(45) COLLATE utf8mb4_unicode_ci NOT NULL COMMENT 'The type of the file.',
  PRIMARY KEY (`idFileTypes`,`type`),
  UNIQUE KEY `idFileTypes_UNIQUE` (`idFileTypes`),
  UNIQUE KEY `FileTypes_type_UNIQUE` (`type`)
) ENGINE=InnoDB AUTO_INCREMENT=10 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `FileTypes`
--

LOCK TABLES `FileTypes` WRITE;
/*!40000 ALTER TABLE `FileTypes` DISABLE KEYS */;
INSERT INTO `FileTypes` VALUES (1,'EXOME_READS'),(2,'GENOME_READS'),(3,'SOMATIC_VCF'),(4,'GERMLINE_VCF'),(5,'EXOME_READS_INDEX'),(6,'GENOME_READS_INDEX'),(7,'SOMATIC_VCF_INDEX'),(8,'GERMLINE_VCF_INDEX'),(9,'UNSTRUCTURED_TEXT');
/*!40000 ALTER TABLE `FileTypes` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `InputFiles`
--

DROP TABLE IF EXISTS `InputFiles`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `InputFiles` (
  `idInputFiles` int(11) unsigned NOT NULL AUTO_INCREMENT COMMENT 'Unique ID of each file.',
  `path` varchar(2048) COLLATE utf8mb4_unicode_ci NOT NULL COMMENT 'Path to the file.',
  `type` varchar(45) COLLATE utf8mb4_unicode_ci NOT NULL COMMENT 'Type of each file.',
  `md5sum` char(32) COLLATE utf8mb4_unicode_ci NOT NULL,
  PRIMARY KEY (`idInputFiles`),
  UNIQUE KEY `idInputFiles_UNIQUE` (`idInputFiles`),
  UNIQUE KEY `InputFiles_md5sum_UNIQUE` (`md5sum`),
  KEY `InputFiles_fileType_idx` (`type`),
  CONSTRAINT `InputFiles_fileType` FOREIGN KEY (`type`) REFERENCES `FileTypes` (`type`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=42 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `InputFiles`
--

LOCK TABLES `InputFiles` WRITE;
/*!40000 ALTER TABLE `InputFiles` DISABLE KEYS */;
INSERT INTO `InputFiles` VALUES (1,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-359781.bai','EXOME_READS_INDEX','66fb4afaeb9bcf1fb6ba0c47e302925b'),(2,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-359781.bam','EXOME_READS','2736eb455d31ede2591aae8216fb6eda'),(3,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-359781.vcf.gz','SOMATIC_VCF','dcec297c7bcd6f3f84566960a27c8b39'),(4,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-359781.vcf.gz.tbi','SOMATIC_VCF_INDEX','d78c8b92abb674356f24e20949dfeba8'),(5,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-359877.bai','EXOME_READS_INDEX','19109d326717b2f16ba4f426483b9a05'),(6,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-359877.bam','EXOME_READS','c9308aed697ef83cb8bfae623cf129f4'),(7,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-359877.vcf.gz','SOMATIC_VCF','d3e640c0c89a34ac7c984ea1772a5b9d'),(8,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-359877.vcf.gz.tbi','SOMATIC_VCF_INDEX','40dc97fe61be9d63b6b7cd9634e451bf'),(9,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-360361.bai','EXOME_READS_INDEX','0ac2525ec8f721e41a07517531b43162'),(10,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-360361.bam','EXOME_READS','34e3ac44bbf74b4012b5c9dc6ac1870d'),(11,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-360361.vcf.gz','SOMATIC_VCF','9c048ce896f1d80923cc1513ce29b137'),(12,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-360361.vcf.gz.tbi','SOMATIC_VCF_INDEX','5e2190bf1c9bdbf8ffafcf31da99ea20'),(13,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-360457.bai','EXOME_READS_INDEX','36e3ff23b02f2e308ebf90dcdeb902d9'),(14,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-360457.bam','EXOME_READS','3646f9c52ab317ee2bc909e28f273f9b'),(15,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-360457.vcf.gz','SOMATIC_VCF','5882c66818f7d5fe09104c485621775b'),(16,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-360457.vcf.gz.tbi','SOMATIC_VCF_INDEX','221fdceeebfa5edaefac873d8cc1905c'),(17,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-361337.bai','EXOME_READS_INDEX','5d8d45855f73b8d066706fa1f211b2b6'),(18,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-361337.bam','EXOME_READS','077d370374b2821defc1fb5d3f85966d'),(19,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-361337.vcf.gz','SOMATIC_VCF','7e5f03272f76f9b18f64b6609506bc76'),(20,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-361337.vcf.gz.tbi','SOMATIC_VCF_INDEX','395408a6c770d20255cb5b11d4ab35e7'),(21,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-361433.bai','EXOME_READS_INDEX','20f2438239dffa9673d35dff58b1381d'),(22,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-361433.bam','EXOME_READS','cc022a43a5454ab956a65a164bc047e4'),(23,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-361433.vcf.gz','SOMATIC_VCF','658583c40648a89a527fe161f5341e7b'),(24,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-361433.vcf.gz.tbi','SOMATIC_VCF_INDEX','ec7a13fea727be1636d38e1af8149d94'),(25,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-362428.bai','EXOME_READS_INDEX','53e700aafd1bc7bddd9a7bc4d4f937e7'),(26,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-362428.bam','EXOME_READS','5a6efabedc3a808ef487429031c8786c'),(27,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-362428.vcf.gz','SOMATIC_VCF','8f0afe2ec90d038271e41d543b34ee7a'),(28,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-362428.vcf.gz.tbi','SOMATIC_VCF_INDEX','40649c2f4baba87710901f1f112fc985'),(29,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-363907.bai','EXOME_READS_INDEX','70578295ea56edd1ce2a62ba41959e6d'),(30,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-363907.bam','EXOME_READS','bb4c81559e3002fcdd09ac48325aecbb'),(31,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-363907.vcf.gz','SOMATIC_VCF','c93f8294c25d1965a810a16bd500be95'),(32,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-363907.vcf.gz.tbi','SOMATIC_VCF_INDEX','7edba1f10be111e128ac726c04dd2ff7'),(33,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-445394.bai','EXOME_READS_INDEX','e4e3336d659b215a91cb81a1a18f4bf9'),(34,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-445394.bam','EXOME_READS','583aa5f341992912221ed3091ff1daf6'),(35,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-445394.vcf.gz','SOMATIC_VCF','880223730fa30dcb60801440bb6a226a'),(36,'gs://broad-dsp-methods-regression-testing/inputData/NexPond-445394.vcf.gz.tbi','SOMATIC_VCF_INDEX','0440125418ac9e3f601b65ad674f0727'),(37,'gs://broad-dsp-methods-regression-testing/inputData/G96830.NA12878.bam','GENOME_READS','7bb97f2a2c6ac23a52ab71d0175d897e'),(38,'gs://broad-dsp-methods-regression-testing/inputData/G96830.NA12878.bai','GENOME_READS_INDEX','9cb49ac1335443baa381e777dd9f485a'),(39,'gs://broad-dsp-methods-regression-testing/inputData/G96830.NA12878.vcf.gz','SOMATIC_VCF','ee32e9a04311032539e16c3ad41a4ea1'),(40,'gs://broad-dsp-methods-regression-testing/inputData/G96830.NA12878.vcf.gz.tbi','SOMATIC_VCF_INDEX','f3bb1d7aa99b9dfa7d82cfd7c3cf39be');
/*!40000 ALTER TABLE `InputFiles` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `Metric_Concordance_FilterAnalysis`
--

DROP TABLE IF EXISTS `Metric_Concordance_FilterAnalysis`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Metric_Concordance_FilterAnalysis` (
  `idConcordance_FilterAnalysis` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `filter` int(11) DEFAULT NULL,
  `tn` int(11) DEFAULT NULL,
  `fn` int(11) DEFAULT NULL,
  `uniqueTn` int(11) DEFAULT NULL,
  `uniqueFn` int(11) DEFAULT NULL,
  PRIMARY KEY (`idConcordance_FilterAnalysis`),
  UNIQUE KEY `idConcordance_FilterAnalysis_UNIQUE` (`idConcordance_FilterAnalysis`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Metric_Concordance_FilterAnalysis`
--

LOCK TABLES `Metric_Concordance_FilterAnalysis` WRITE;
/*!40000 ALTER TABLE `Metric_Concordance_FilterAnalysis` DISABLE KEYS */;
/*!40000 ALTER TABLE `Metric_Concordance_FilterAnalysis` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `Metric_Concordance_Summary`
--

DROP TABLE IF EXISTS `Metric_Concordance_Summary`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Metric_Concordance_Summary` (
  `idConcordance_Summary` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `type` varchar(512) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `truePositive` double DEFAULT NULL,
  `falsePositive` double DEFAULT NULL,
  `falseNegative` double DEFAULT NULL,
  `sensitivity` double DEFAULT NULL,
  `precision` double DEFAULT NULL,
  PRIMARY KEY (`idConcordance_Summary`),
  UNIQUE KEY `idConcordance_Summary_UNIQUE` (`idConcordance_Summary`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Metric_Concordance_Summary`
--

LOCK TABLES `Metric_Concordance_Summary` WRITE;
/*!40000 ALTER TABLE `Metric_Concordance_Summary` DISABLE KEYS */;
/*!40000 ALTER TABLE `Metric_Concordance_Summary` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `Metric_GenotypeConcordance_ContingencyMetrics`
--

DROP TABLE IF EXISTS `Metric_GenotypeConcordance_ContingencyMetrics`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Metric_GenotypeConcordance_ContingencyMetrics` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `variantType` varchar(256) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `truthSample` varchar(256) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `callSample` varchar(256) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `tpCount` int(11) DEFAULT NULL,
  `tnCount` int(11) DEFAULT NULL,
  `fpCount` int(11) DEFAULT NULL,
  `fnCount` int(11) DEFAULT NULL,
  `emptyCount` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `id_metric_GenotypeConcordance_ContingencyMetrics_UNIQUE` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Metric_GenotypeConcordance_ContingencyMetrics`
--

LOCK TABLES `Metric_GenotypeConcordance_ContingencyMetrics` WRITE;
/*!40000 ALTER TABLE `Metric_GenotypeConcordance_ContingencyMetrics` DISABLE KEYS */;
/*!40000 ALTER TABLE `Metric_GenotypeConcordance_ContingencyMetrics` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `Metric_GenotypeConcordance_DetailMetrics`
--

DROP TABLE IF EXISTS `Metric_GenotypeConcordance_DetailMetrics`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Metric_GenotypeConcordance_DetailMetrics` (
  `idGenotypeConcordance_DetailMetrics` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `variantType` varchar(512) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `truthSample` varchar(512) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `callSample` varchar(512) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `truthState` varchar(512) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `callState` varchar(512) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `count` int(11) DEFAULT NULL,
  `contingencyValues` varchar(512) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`idGenotypeConcordance_DetailMetrics`),
  UNIQUE KEY `idGenotypeConcordance_DetailMetrics_UNIQUE` (`idGenotypeConcordance_DetailMetrics`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Metric_GenotypeConcordance_DetailMetrics`
--

LOCK TABLES `Metric_GenotypeConcordance_DetailMetrics` WRITE;
/*!40000 ALTER TABLE `Metric_GenotypeConcordance_DetailMetrics` DISABLE KEYS */;
/*!40000 ALTER TABLE `Metric_GenotypeConcordance_DetailMetrics` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `Metric_GenotypeConcordance_SummaryMetrics`
--

DROP TABLE IF EXISTS `Metric_GenotypeConcordance_SummaryMetrics`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Metric_GenotypeConcordance_SummaryMetrics` (
  `idGenotypeConcordance_SummaryMetrics` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `variantType` varchar(512) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `truthSample` varchar(512) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `callSample` varchar(512) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `hetSensitivity` double DEFAULT NULL,
  `hetPPV` double DEFAULT NULL,
  `hetSpecificity` double DEFAULT NULL,
  `homvarSensitivity` double DEFAULT NULL,
  `varSensitivity` double DEFAULT NULL,
  `varPPV` double DEFAULT NULL,
  `varSpecificity` double DEFAULT NULL,
  `genotypeConcordance` double DEFAULT NULL,
  `nonRefGenotypeConcordance` double DEFAULT NULL,
  PRIMARY KEY (`idGenotypeConcordance_SummaryMetrics`),
  UNIQUE KEY `idGenotypeConcordance_SummaryMetrics_UNIQUE` (`idGenotypeConcordance_SummaryMetrics`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Metric_GenotypeConcordance_SummaryMetrics`
--

LOCK TABLES `Metric_GenotypeConcordance_SummaryMetrics` WRITE;
/*!40000 ALTER TABLE `Metric_GenotypeConcordance_SummaryMetrics` DISABLE KEYS */;
/*!40000 ALTER TABLE `Metric_GenotypeConcordance_SummaryMetrics` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `Metric_Timing`
--

DROP TABLE IF EXISTS `Metric_Timing`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Metric_Timing` (
  `idMetric_Timing` int(11) unsigned NOT NULL,
  `startTime` double DEFAULT NULL COMMENT 'Epoch start time of this tool.',
  `endTime` double DEFAULT NULL COMMENT 'Epoch end time of this tool.',
  `elapsedTime` double DEFAULT NULL COMMENT 'Elapsed time of this tool (seconds).',
  PRIMARY KEY (`idMetric_Timing`),
  UNIQUE KEY `idMetric_Timing_UNIQUE` (`idMetric_Timing`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Metric_Timing`
--

LOCK TABLES `Metric_Timing` WRITE;
/*!40000 ALTER TABLE `Metric_Timing` DISABLE KEYS */;
/*!40000 ALTER TABLE `Metric_Timing` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `Metrics`
--

DROP TABLE IF EXISTS `Metrics`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Metrics` (
  `idMetrics` int(11) unsigned NOT NULL,
  `metricTableName` varchar(1024) COLLATE utf8mb4_unicode_ci NOT NULL COMMENT 'The name of the metric table containing the information.  Also the name of the metric collected.',
  `concreteMetricID` int(11) NOT NULL,
  `sourceAnalysis` int(11) unsigned NOT NULL,
  PRIMARY KEY (`idMetrics`),
  UNIQUE KEY `idMetrics_UNIQUE` (`idMetrics`),
  KEY `Metrics_sourceAnalysisID` (`sourceAnalysis`),
  CONSTRAINT `Metrics_sourceAnalysisID` FOREIGN KEY (`sourceAnalysis`) REFERENCES `AnalysisRuns` (`idAnalysisRun`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Metrics`
--

LOCK TABLES `Metrics` WRITE;
/*!40000 ALTER TABLE `Metrics` DISABLE KEYS */;
/*!40000 ALTER TABLE `Metrics` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `MetricsEvaluator`
--

DROP TABLE IF EXISTS `MetricsEvaluator`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `MetricsEvaluator` (
  `idMetricsEvaluator` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `scenarioInfo` int(11) unsigned NOT NULL,
  `metric` int(11) unsigned NOT NULL,
  `absoluteAllowableMaxValue` double DEFAULT NULL,
  `absoluteAllowableMinValue` double DEFAULT NULL,
  `pairwiseAllowableIncrease` double DEFAULT NULL,
  `pairwiseAllowableDecrease` double DEFAULT NULL,
  `isPairwiseMeasurePercentage` tinyint(4) DEFAULT NULL,
  PRIMARY KEY (`idMetricsEvaluator`),
  UNIQUE KEY `idMetricsEvaluator_UNIQUE` (`idMetricsEvaluator`),
  KEY `MetricsEvaluator_metricID_idx` (`metric`),
  KEY `MetricsEvaluator_scenarioInfoID_idx` (`scenarioInfo`),
  CONSTRAINT `MetricsEvaluator_metricID` FOREIGN KEY (`metric`) REFERENCES `Metrics` (`idMetrics`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `MetricsEvaluator_scenarioInfoID` FOREIGN KEY (`scenarioInfo`) REFERENCES `TestScenarioInfo` (`idTestScenarioInfo`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `MetricsEvaluator`
--

LOCK TABLES `MetricsEvaluator` WRITE;
/*!40000 ALTER TABLE `MetricsEvaluator` DISABLE KEYS */;
/*!40000 ALTER TABLE `MetricsEvaluator` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `MetricsEvaluatorData`
--

DROP TABLE IF EXISTS `MetricsEvaluatorData`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `MetricsEvaluatorData` (
  `metricsEvaluatorID` int(11) unsigned NOT NULL,
  `dataID` int(11) unsigned NOT NULL,
  PRIMARY KEY (`metricsEvaluatorID`,`dataID`),
  KEY `MetricsEvaluatorData_dataID_idx` (`dataID`),
  CONSTRAINT `MetricsEvaluatorData_dataID` FOREIGN KEY (`dataID`) REFERENCES `OutputFiles` (`idOutputFiles`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `MetricsEvaluatorData_metricsEvaluatorID` FOREIGN KEY (`metricsEvaluatorID`) REFERENCES `MetricsEvaluator` (`idMetricsEvaluator`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `MetricsEvaluatorData`
--

LOCK TABLES `MetricsEvaluatorData` WRITE;
/*!40000 ALTER TABLE `MetricsEvaluatorData` DISABLE KEYS */;
/*!40000 ALTER TABLE `MetricsEvaluatorData` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `OutputFiles`
--

DROP TABLE IF EXISTS `OutputFiles`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `OutputFiles` (
  `idOutputFiles` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `fileType` varchar(45) COLLATE utf8mb4_unicode_ci NOT NULL,
  `path` varchar(2048) COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `timeCreated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `md5sum` char(32) COLLATE utf8mb4_unicode_ci NOT NULL,
  `sourceType` enum('TOOL','ANALYSIS') COLLATE utf8mb4_unicode_ci NOT NULL,
  PRIMARY KEY (`idOutputFiles`),
  UNIQUE KEY `idOutputFiles_UNIQUE` (`idOutputFiles`),
  KEY `OutputFiles_fileType_idx` (`fileType`),
  CONSTRAINT `OutputFiles_fileType` FOREIGN KEY (`fileType`) REFERENCES `FileTypes` (`type`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `OutputFiles`
--

LOCK TABLES `OutputFiles` WRITE;
/*!40000 ALTER TABLE `OutputFiles` DISABLE KEYS */;
/*!40000 ALTER TABLE `OutputFiles` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `TestReport`
--

DROP TABLE IF EXISTS `TestReport`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `TestReport` (
  `idTestReport` int(11) unsigned NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`idTestReport`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `TestReport`
--

LOCK TABLES `TestReport` WRITE;
/*!40000 ALTER TABLE `TestReport` DISABLE KEYS */;
/*!40000 ALTER TABLE `TestReport` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `TestReportData`
--

DROP TABLE IF EXISTS `TestReportData`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `TestReportData` (
  `reportID` int(11) unsigned NOT NULL,
  `dataID` int(11) unsigned NOT NULL,
  PRIMARY KEY (`reportID`,`dataID`),
  KEY `TestReportData_dataID_idx` (`dataID`),
  CONSTRAINT `TestReportData_dataID` FOREIGN KEY (`dataID`) REFERENCES `OutputFiles` (`idOutputFiles`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `TestReportData_reportID` FOREIGN KEY (`reportID`) REFERENCES `TestReport` (`idTestReport`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `TestReportData`
--

LOCK TABLES `TestReportData` WRITE;
/*!40000 ALTER TABLE `TestReportData` DISABLE KEYS */;
/*!40000 ALTER TABLE `TestReportData` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `TestReportMetrics`
--

DROP TABLE IF EXISTS `TestReportMetrics`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `TestReportMetrics` (
  `reportId` int(11) unsigned NOT NULL,
  `metricId` int(11) unsigned NOT NULL,
  KEY `TestReportMetrics_reportID_idx` (`reportId`),
  KEY `TestReportMetrics_metricID_idx` (`metricId`),
  CONSTRAINT `TestReportMetrics_metricID` FOREIGN KEY (`metricId`) REFERENCES `Metrics` (`idMetrics`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `TestReportMetrics_reportID` FOREIGN KEY (`reportId`) REFERENCES `TestReport` (`idTestReport`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `TestReportMetrics`
--

LOCK TABLES `TestReportMetrics` WRITE;
/*!40000 ALTER TABLE `TestReportMetrics` DISABLE KEYS */;
/*!40000 ALTER TABLE `TestReportMetrics` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `TestScenario`
--

DROP TABLE IF EXISTS `TestScenario`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `TestScenario` (
  `idTestScenarioRun` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `scenarioInfo` int(11) unsigned NOT NULL,
  `configuration` json NOT NULL,
  `inputFile` int(11) unsigned NOT NULL,
  `tool` int(11) unsigned NOT NULL,
  PRIMARY KEY (`idTestScenarioRun`),
  UNIQUE KEY `idTestScenario_UNIQUE` (`idTestScenarioRun`),
  KEY `TestScenario_fileID_idx` (`inputFile`),
  KEY `TestScenario_toolID_idx` (`tool`),
  KEY `TestScenario_scenarioInfoID_idx` (`scenarioInfo`),
  CONSTRAINT `TestScenario_fileID` FOREIGN KEY (`inputFile`) REFERENCES `InputFiles` (`idInputFiles`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `TestScenario_scenarioInfoID` FOREIGN KEY (`scenarioInfo`) REFERENCES `TestScenarioInfo` (`idTestScenarioInfo`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `TestScenario_toolID` FOREIGN KEY (`tool`) REFERENCES `Tools` (`idTools`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=100 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `TestScenario`
--

LOCK TABLES `TestScenario` WRITE;
/*!40000 ALTER TABLE `TestScenario` DISABLE KEYS */;
INSERT INTO `TestScenario` VALUES (1,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359781.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-07-4.1.0.0-50-g342569ac0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',2,1),(2,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359877.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-07-4.1.0.0-50-g342569ac0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',6,1),(3,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360361.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-07-4.1.0.0-50-g342569ac0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',10,1),(4,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360457.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-07-4.1.0.0-50-g342569ac0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',14,1),(5,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361337.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-07-4.1.0.0-50-g342569ac0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',18,1),(6,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361433.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-07-4.1.0.0-50-g342569ac0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',22,1),(7,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-362428.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-07-4.1.0.0-50-g342569ac0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',26,1),(8,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-363907.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-07-4.1.0.0-50-g342569ac0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',30,1),(9,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-445394.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-07-4.1.0.0-50-g342569ac0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',34,1),(10,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359781.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.4.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',2,2),(11,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359877.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.4.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',6,2),(12,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360361.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.4.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',10,2),(13,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360457.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.4.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',14,2),(14,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361337.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.4.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',18,2),(15,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361433.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.4.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',22,2),(16,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-362428.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.4.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',26,2),(17,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-363907.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.4.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',30,2),(18,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-445394.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.4.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',34,2),(19,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359781.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.6.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',2,3),(20,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359877.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.6.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',6,3),(21,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360361.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.6.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',10,3),(22,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360457.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.6.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',14,3),(23,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361337.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.6.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',18,3),(24,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361433.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.6.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',22,3),(25,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-362428.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.6.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',26,3),(26,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-363907.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.6.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',30,3),(27,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-445394.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.6.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',34,3),(28,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359781.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.7.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',2,4),(29,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359877.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.7.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',6,4),(30,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360361.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.7.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',10,4),(31,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360457.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.7.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',14,4),(32,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361337.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.7.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',18,4),(33,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361433.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.7.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',22,4),(34,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-362428.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.7.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',26,4),(35,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-363907.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.7.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',30,4),(36,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-445394.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.7.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',34,4),(37,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359781.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',2,5),(38,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359877.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',6,5),(39,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360361.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',10,5),(40,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360457.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',14,5),(41,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361337.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',18,5),(42,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361433.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',22,5),(43,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-362428.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',26,5),(44,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-363907.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',30,5),(45,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-445394.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',34,5),(46,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359781.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',2,6),(47,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359877.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',6,6),(48,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360361.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',10,6),(49,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360457.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',14,6),(50,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361337.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',18,6),(51,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361433.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',22,6),(52,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-362428.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',26,6),(53,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-363907.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',30,6),(54,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-445394.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.8.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',34,6),(55,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359781.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.9.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',2,7),(56,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359877.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.9.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',6,7),(57,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360361.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.9.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',10,7),(58,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360457.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.9.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',14,7),(59,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361337.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.9.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',18,7),(60,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361433.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.9.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',22,7),(61,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-362428.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.9.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',26,7),(62,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-363907.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.9.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',30,7),(63,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-445394.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.9.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',34,7),(64,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359781.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',2,8),(65,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359877.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',6,8),(66,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360361.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',10,8),(67,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360457.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',14,8),(68,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361337.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',18,8),(69,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361433.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',22,8),(70,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-362428.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',26,8),(71,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-363907.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',30,8),(72,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-445394.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.1-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',34,8),(73,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359781.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',2,9),(74,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359877.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',6,9),(75,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360361.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',10,9),(76,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360457.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',14,9),(77,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361337.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',18,9),(78,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361433.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',22,9),(79,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-362428.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',26,9),(80,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-363907.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',30,9),(81,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-445394.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.10.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',34,9),(82,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359781.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.12.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',2,10),(83,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359877.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.12.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',6,10),(84,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360361.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.12.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',10,10),(85,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360457.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.12.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',14,10),(86,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361337.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.12.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',18,10),(87,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361433.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.12.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',22,10),(88,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-362428.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.12.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',26,10),(89,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-363907.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.12.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',30,10),(90,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-445394.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.0.12.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',34,10),(91,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359781.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.1.0.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',2,11),(92,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-359877.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.1.0.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',6,11),(93,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360361.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.1.0.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',10,11),(94,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-360457.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.1.0.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',14,11),(95,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361337.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.1.0.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',18,11),(96,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-361433.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.1.0.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',22,11),(97,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-362428.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.1.0.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',26,11),(98,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-363907.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.1.0.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',30,11),(99,1,'{\"HaplotypeCallerTask.mem_gb\": 32, \"HaplotypeCallerTask.gvcf_mode\": \"false\", \"HaplotypeCallerTask.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\", \"HaplotypeCallerTask.input_bams\": [\"NexPond-445394.bam\"], \"HaplotypeCallerTask.gatk_docker\": \"jonnsmith/gatk_test_builds:2019-03-08-4.1.0.0-SNAPSHOT\", \"HaplotypeCallerTask.contamination\": 0, \"HaplotypeCallerTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\", \"HaplotypeCallerTask.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\", \"HaplotypeCallerTask.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\", \"HaplotypeCallerTask.interval_padding\": 0, \"HaplotypeCallerTask.boot_disk_size_gb\": 64, \"HaplotypeCallerTask.default_disk_space_gb\": 512}',34,11);
/*!40000 ALTER TABLE `TestScenario` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `TestScenarioInfo`
--

DROP TABLE IF EXISTS `TestScenarioInfo`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `TestScenarioInfo` (
  `idTestScenarioInfo` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(512) COLLATE utf8mb4_unicode_ci NOT NULL,
  `description` varchar(2048) COLLATE utf8mb4_unicode_ci NOT NULL,
  PRIMARY KEY (`idTestScenarioInfo`),
  UNIQUE KEY `idTestScenarioInfo_UNIQUE` (`idTestScenarioInfo`),
  UNIQUE KEY `name_UNIQUE` (`name`)
) ENGINE=InnoDB AUTO_INCREMENT=2 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `TestScenarioInfo`
--

LOCK TABLES `TestScenarioInfo` WRITE;
/*!40000 ALTER TABLE `TestScenarioInfo` DISABLE KEYS */;
INSERT INTO `TestScenarioInfo` VALUES (1,'HaplotypeCaller Validation','Validation tests for the HaplotypeCaller.  Includes both correctness and performance metrics.');
/*!40000 ALTER TABLE `TestScenarioInfo` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `TestScenarioOutputFiles`
--

DROP TABLE IF EXISTS `TestScenarioOutputFiles`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `TestScenarioOutputFiles` (
  `scenarioRun` int(11) unsigned NOT NULL,
  `outputFile` int(11) unsigned NOT NULL,
  PRIMARY KEY (`scenarioRun`,`outputFile`),
  KEY `TestScenarioOutputFiles_outputFileID_idx` (`outputFile`),
  CONSTRAINT `TestScenarioOutputFiles_outputFileID` FOREIGN KEY (`outputFile`) REFERENCES `OutputFiles` (`idOutputFiles`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `TestScenarioOutputFiles_scenarioRunID` FOREIGN KEY (`scenarioRun`) REFERENCES `TestScenario` (`idTestScenarioRun`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `TestScenarioOutputFiles`
--

LOCK TABLES `TestScenarioOutputFiles` WRITE;
/*!40000 ALTER TABLE `TestScenarioOutputFiles` DISABLE KEYS */;
/*!40000 ALTER TABLE `TestScenarioOutputFiles` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `Tools`
--

DROP TABLE IF EXISTS `Tools`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Tools` (
  `idTools` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `version` varchar(1024) COLLATE utf8mb4_unicode_ci NOT NULL,
  `name` varchar(1024) COLLATE utf8mb4_unicode_ci NOT NULL,
  `dateReleased` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `wdl` varchar(1024) COLLATE utf8mb4_unicode_ci NOT NULL,
  `wdlChecksum` char(32) COLLATE utf8mb4_unicode_ci NOT NULL,
  PRIMARY KEY (`idTools`),
  UNIQUE KEY `idTools_UNIQUE` (`idTools`)
) ENGINE=InnoDB AUTO_INCREMENT=12 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Tools`
--

LOCK TABLES `Tools` WRITE;
/*!40000 ALTER TABLE `Tools` DISABLE KEYS */;
INSERT INTO `Tools` VALUES (1,'2019-03-07-4.1.0.0-50-g342569ac0-SNAPSHOT','HaplotypeCaller','2019-03-07 21:21:30','gs://broad-dsp-methods-regression-testing/toolWdls/haplotypeCaller.wdl','30c737cea0697e08eeb79214a068fa97'),(2,'2019-03-08-4.0.4.0-SNAPSHOT','HaplotypeCaller','2018-04-26 13:46:04','gs://broad-dsp-methods-regression-testing/toolWdls/haplotypeCaller.wdl','30c737cea0697e08eeb79214a068fa97'),(3,'2019-03-08-4.0.6.0-SNAPSHOT','HaplotypeCaller','2018-07-06 21:40:31','gs://broad-dsp-methods-regression-testing/toolWdls/haplotypeCaller.wdl','30c737cea0697e08eeb79214a068fa97'),(4,'2019-03-08-4.0.7.0-SNAPSHOT','HaplotypeCaller','2018-07-30 13:46:28','gs://broad-dsp-methods-regression-testing/toolWdls/haplotypeCaller.wdl','30c737cea0697e08eeb79214a068fa97'),(5,'2019-03-08-4.0.8.0-SNAPSHOT','HaplotypeCaller','2018-08-14 18:50:47','gs://broad-dsp-methods-regression-testing/toolWdls/haplotypeCaller.wdl','30c737cea0697e08eeb79214a068fa97'),(6,'2019-03-08-4.0.8.1-SNAPSHOT','HaplotypeCaller','2018-08-16 19:03:52','gs://broad-dsp-methods-regression-testing/toolWdls/haplotypeCaller.wdl','30c737cea0697e08eeb79214a068fa97'),(7,'2019-03-08-4.0.9.0-SNAPSHOT','HaplotypeCaller','2018-09-20 02:27:12','gs://broad-dsp-methods-regression-testing/toolWdls/haplotypeCaller.wdl','30c737cea0697e08eeb79214a068fa97'),(8,'2019-03-08-4.0.10.1-SNAPSHOT','HaplotypeCaller','2018-10-03 20:11:05','gs://broad-dsp-methods-regression-testing/toolWdls/haplotypeCaller.wdl','30c737cea0697e08eeb79214a068fa97'),(9,'2019-03-08-4.0.10.0-SNAPSHOT','HaplotypeCaller','2018-10-09 17:33:58','gs://broad-dsp-methods-regression-testing/toolWdls/haplotypeCaller.wdl','30c737cea0697e08eeb79214a068fa97'),(10,'2019-03-08-4.0.12.0-SNAPSHOT','HaplotypeCaller','2018-12-17 14:47:45','gs://broad-dsp-methods-regression-testing/toolWdls/haplotypeCaller.wdl','30c737cea0697e08eeb79214a068fa97'),(11,'2019-03-08-4.1.0.0-SNAPSHOT','HaplotypeCaller','2019-01-30 02:44:18','gs://broad-dsp-methods-regression-testing/toolWdls/haplotypeCaller.wdl','30c737cea0697e08eeb79214a068fa97');
/*!40000 ALTER TABLE `Tools` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2019-03-14 14:30:31
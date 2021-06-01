
pipeline {
  agent any
  stages {
    stage('build') {
      when {
        not {
          buildingTag()
        }
      }
      steps {
        sh '/.mounts/labs/gsi/vidarr/jenkins-ci-wrapper test -t /.mounts/labs/gsi/vidarr/testing-config.json'
      }
    }
    stage('Deploy') {
      when {
        buildingTag()
      }
      steps {
        sh '/.mounts/labs/gsi/vidarr/jenkins-ci-wrapper deploy -v $TAG_NAME -t /.mounts/labs/gsi/vidarr/testing-config.json -U /.mounts/labs/gsi/vidarr/deploy-urls'
      }
    }
  }
}

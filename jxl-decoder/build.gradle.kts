plugins {
    kotlin("jvm")
    `java-library`
    `maven-publish`
}

kotlin {
    // Android 12 (API 31) 以降の ART でも動くように Java 17 のバイトコードにする
    compilerOptions {
        jvmTarget.set(org.jetbrains.kotlin.gradle.dsl.JvmTarget.JVM_17)
        // CI の lint では -Pwerror を付けて警告をエラーにする
        allWarningsAsErrors.set(providers.gradleProperty("werror").isPresent)
    }
    explicitApi()
}

java {
    sourceCompatibility = JavaVersion.VERSION_17
    targetCompatibility = JavaVersion.VERSION_17
    withSourcesJar()
}

dependencies {
    testImplementation(kotlin("test"))
}

// ライセンスを jar (クラスとソースの両方) に入れる
tasks.withType<Jar>().configureEach {
    from(rootProject.files("LICENSE", "LICENSE-jxl-oxide")) { into("META-INF") }
}

tasks.test {
    useJUnitPlatform()
    maxHeapSize = "2g"
    systemProperty("jxl.testdata", rootProject.file("testdata").absolutePath)
}

publishing {
    publications {
        create<MavenPublication>("maven") {
            artifactId = "jxl-decoder"
            from(components["java"])
        }
    }
}


// GitHub のリリースを Ivy リポジトリとして使うための ivy-<version>.xml (CI がリリースに添付する)。
// jar だけだと Gradle は HEAD で存在を確かめるが、GitHub のリリースの転送先は HEAD で失敗することがあるので、
// まずこのファイルを GET で読ませる
val ivyDescriptor by tasks.registering {
    val version = project.version.toString()
    val kotlinVersion = org.jetbrains.kotlin.gradle.plugin.getKotlinPluginVersion(logger)
    val out = layout.buildDirectory.file("ivy/ivy-$version.xml")
    inputs.property("version", version)
    inputs.property("kotlinVersion", kotlinVersion)
    outputs.file(out)
    doLast {
        out.get().asFile.writeText(
            """
            <?xml version="1.0" encoding="UTF-8"?>
            <ivy-module version="2.0" xmlns:m="http://ant.apache.org/ivy/maven">
              <info organisation="io.github.hamtit" module="jxl-decoder" revision="$version" status="release"/>
              <configurations>
                <conf name="default"/>
                <conf name="sources"/>
              </configurations>
              <publications>
                <artifact name="jxl-decoder" type="jar" ext="jar" conf="default"/>
                <artifact name="jxl-decoder" type="source" ext="jar" conf="sources" m:classifier="sources"/>
              </publications>
              <dependencies>
                <dependency org="org.jetbrains.kotlin" name="kotlin-stdlib" rev="$kotlinVersion" conf="default->default"/>
              </dependencies>
            </ivy-module>
            """.trimIndent() + "\n",
        )
    }
}
